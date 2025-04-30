#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import re
import gzip
import os
import subprocess
import tempfile
import requests
import io
import shutil
from functools import partial
from collections import defaultdict
from cyvcf2 import VCF, Writer
from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl import Workbook
import glob

def download_feature_table(reference_id, output_dir="tables", show_progress=False):
    """
    Скачивает таблицу признаков для указанного референса из NCBI.
    
    Функция использует структуру FTP сервера NCBI для поиска и загрузки
    таблицы признаков генома.
    """
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, "feature_table.tsv")
    
    # Метод 1: Поиск через E-utilities API для получения правильного accession номера
    try:
        if show_progress:
            print(f"Поиск информации о геноме {reference_id} через NCBI...")
            
        # Подготовка поискового запроса
        if "_" in reference_id:
            parts = reference_id.split("_")
            genus = parts[0]
            if len(parts) > 1:
                species = "_".join(parts[1:])
                search_term = f"{genus}+{species}"
            else:
                search_term = reference_id
        else:
            search_term = reference_id
            
        # Поиск в базе данных Assembly
        search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=assembly&term={search_term}&retmax=10&retmode=json"
        search_response = requests.get(search_url)
        search_data = search_response.json()
        
        if "esearchresult" in search_data and "idlist" in search_data["esearchresult"] and search_data["esearchresult"]["idlist"]:
            assembly_ids = search_data["esearchresult"]["idlist"]
            
            # Получаем подробную информацию о найденных сборках
            for assembly_id in assembly_ids:
                assembly_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={assembly_id}&retmode=json"
                assembly_response = requests.get(assembly_url)
                assembly_data = assembly_response.json()
                
                if "result" in assembly_data and assembly_id in assembly_data["result"]:
                    assembly_info = assembly_data["result"][assembly_id]
                    
                    # Получаем accession номер
                    if "assemblyaccession" in assembly_info:
                        accession = assembly_info["assemblyaccession"]
                        assembly_name = assembly_info.get("assemblyname", "")
                        
                        if show_progress:
                            print(f"Найден геном: {accession} ({assembly_name})")
                        
                        # Формируем URL для FTP доступа по шаблону из примера
                        # Например: GCF_000195955.2 -> GCF/000/195/955/GCF_000195955.2_ASM19595v2
                        
                        # Убираем версию из accession, если она есть
                        base_accession = accession.split('.')[0] if '.' in accession else accession
                        
                        # Разбиваем числовую часть на группы по 3 цифры
                        numeric_part = base_accession.split('_')[1]
                        path_parts = [numeric_part[i:i+3] for i in range(0, len(numeric_part), 3)]
                        
                        # Формируем путь в FTP структуре
                        prefix = base_accession.split('_')[0]  # GCF или GCA
                        
                        # Формируем возможные варианты полного имени каталога
                        possible_dirs = [
                            f"{accession}_{assembly_name.replace(' ', '_')}",
                            f"{accession}_ASM{numeric_part}v2",
                            accession
                        ]
                        
                        # Проверяем каждый вариант
                        for dir_name in possible_dirs:
                            ftp_path = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{'/'.join(path_parts)}/{dir_name}"
                            feature_table_url = f"{ftp_path}/{dir_name}_feature_table.txt.gz"
                            
                            if show_progress:
                                print(f"Проверяем URL: {feature_table_url}")
                            
                            try:
                                response = requests.head(feature_table_url)
                                if response.status_code == 200:
                                    if show_progress:
                                        print(f"Найдена таблица признаков: {feature_table_url}")
                                        print(f"Загрузка...")
                                    
                                    # Загружаем файл и обрабатываем
                                    download_response = requests.get(feature_table_url, stream=True)
                                    download_response.raise_for_status()
                                    
                                    # Сохраняем во временный файл
                                    with tempfile.NamedTemporaryFile(delete=False, suffix='.gz') as temp_file:
                                        for chunk in download_response.iter_content(chunk_size=8192):
                                            if chunk:
                                                temp_file.write(chunk)
                                        temp_file_path = temp_file.name
                                    
                                    # Обрабатываем и преобразуем данные
                                    data = []
                                    with gzip.open(temp_file_path, 'rt') as f:
                                        # Пропускаем заголовок
                                        header = next(f)
                                        
                                        for line in f:
                                            fields = line.strip().split('\t')
                                            if len(fields) >= 15:  # Проверяем наличие достаточного количества полей
                                                feature = fields[0]
                                                locus_tag = fields[6]
                                                symbol = fields[7]
                                                product = fields[14] if len(fields) > 14 else ""
                                                
                                                # Включаем только записи для генов
                                                if feature.lower() == "gene" and locus_tag:
                                                    data.append([symbol, locus_tag, product])
                                    
                                    # Удаляем временный файл
                                    os.unlink(temp_file_path)
                                    
                                    # Сохраняем данные в TSV
                                    if data:
                                        df = pd.DataFrame(data, columns=["symbol", "locus_tag", "name"])
                                        df.to_csv(output_path, sep="\t", index=False)
                                        
                                        if show_progress:
                                            print(f"Таблица признаков успешно загружена ({len(data)} генов) и сохранена в {output_path}")
                                        return output_path
                                    else:
                                        if show_progress:
                                            print(f"Предупреждение: Файл загружен, но не содержит данных о генах.")
                            except Exception as e:
                                print(f"Ошибка при проверке {feature_table_url}: {str(e)}")
                                continue
    
    except Exception as e:
        print(f"Ошибка при поиске через NCBI API: {str(e)}")
    
    # Метод 2: Если не нашли через API, пробуем прямой поиск по известным шаблонам
    try:
        if show_progress:
            print(f"Попытка прямого поиска в каталогах NCBI FTP...")
        
        # Пробуем различные комбинации для accession
        accession_patterns = [
            f"GCF_{reference_id}",
            f"GCA_{reference_id}",
            reference_id
        ]
        
        for acc in accession_patterns:
            # Разбиваем на компоненты
            if "_" in acc:
                prefix, numeric_part = acc.split("_", 1)
                
                # Если есть версия (например, .1 или .2), отделяем её
                if "." in numeric_part:
                    numeric_part = numeric_part.split(".")[0]
                
                # Разбиваем числовую часть на группы по 3 цифры
                path_parts = []
                for i in range(0, len(numeric_part), 3):
                    if i+3 <= len(numeric_part):
                        path_parts.append(numeric_part[i:i+3])
                    else:
                        path_parts.append(numeric_part[i:].ljust(3, '0'))
                
                # Формируем базовый путь на FTP
                ftp_base = f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{prefix}/{'/'.join(path_parts)}"
                
                # Получаем список каталогов в указанном FTP пути
                try:
                    response = requests.get(ftp_base)
                    if response.status_code == 200:
                        # Ищем подходящий каталог
                        href_pattern = re.compile(r'href="([^"]+)"')
                        for match in href_pattern.finditer(response.text):
                            dir_name = match.group(1)
                            if acc in dir_name and dir_name.endswith('/'):
                                # Нашли каталог, проверяем наличие feature_table
                                dir_path = dir_name.rstrip('/')
                                feature_table_url = f"{ftp_base}/{dir_path}/{dir_path}_feature_table.txt.gz"
                                
                                if show_progress:
                                    print(f"Проверяем URL: {feature_table_url}")
                                
                                try:
                                    download_response = requests.head(feature_table_url)
                                    if download_response.status_code == 200:
                                        if show_progress:
                                            print(f"Найдена таблица признаков: {feature_table_url}")
                                            print(f"Загрузка...")
                                        
                                        # Загружаем и обрабатываем файл
                                        feature_response = requests.get(feature_table_url, stream=True)
                                        feature_response.raise_for_status()
                                        
                                        with tempfile.NamedTemporaryFile(delete=False, suffix='.gz') as temp_file:
                                            for chunk in feature_response.iter_content(chunk_size=8192):
                                                if chunk:
                                                    temp_file.write(chunk)
                                            temp_file_path = temp_file.name
                                        
                                        # Обрабатываем данные
                                        data = []
                                        with gzip.open(temp_file_path, 'rt') as f:
                                            # Пропускаем заголовок
                                            header = next(f)
                                            
                                            for line in f:
                                                fields = line.strip().split('\t')
                                                if len(fields) >= 15:
                                                    feature = fields[0]
                                                    locus_tag = fields[6]
                                                    symbol = fields[7]
                                                    product = fields[14] if len(fields) > 14 else ""
                                                    
                                                    if feature.lower() == "gene" and locus_tag:
                                                        data.append([symbol, locus_tag, product])
                                        
                                        # Удаляем временный файл
                                        os.unlink(temp_file_path)
                                        
                                        # Сохраняем данные
                                        if data:
                                            df = pd.DataFrame(data, columns=["symbol", "locus_tag", "name"])
                                            df.to_csv(output_path, sep="\t", index=False)
                                            
                                            if show_progress:
                                                print(f"Таблица признаков успешно загружена ({len(data)} генов) и сохранена в {output_path}")
                                            return output_path
                                except Exception as e:
                                    print(f"Ошибка при проверке {feature_table_url}: {str(e)}")
                                    continue
                except Exception as e:
                    print(f"Ошибка при чтении каталога {ftp_base}: {str(e)}")
                    continue
    
    except Exception as e:
        print(f"Ошибка при прямом поиске в каталогах NCBI FTP: {str(e)}")
    
    # Если все методы не сработали, создаем пустую таблицу
    if show_progress:
        print(f"Не удалось найти или загрузить таблицу признаков. Создаем пустую таблицу.")
    
    df = pd.DataFrame(columns=["symbol", "locus_tag", "name"])
    df.to_csv(output_path, sep="\t", index=False)
    
    print(f"ВНИМАНИЕ: Создана пустая таблица признаков в {output_path}. Для корректной работы скрипта,")
    print(f"вам необходимо заполнить её вручную или загрузить соответствующий файл с аннотациями.")
    print(f"Файл для большинства геномов можно найти вручную на:")
    print(f"https://ftp.ncbi.nlm.nih.gov/genomes/all/")
    print(f"Поиск по идентификатору {reference_id} или соответствующему accession номеру.")
    print(f"Файл должен иметь формат ACCESSION_feature_table.txt.gz")
    print(f"После загрузки, распакуйте его и преобразуйте в требуемый формат:")
    print(f"symbol\tlocus_tag\tname")
    
    return output_path

def ensure_feature_table(reference_id, custom_table=None, show_progress=False):
    """
    Проверяет наличие таблицы признаков и скачивает её при необходимости.
    
    Позволяет использовать пользовательскую таблицу или скачать соответствующую
    таблицу для указанного референса.
    """
    if custom_table and os.path.exists(custom_table):
        if show_progress:
            print(f"Используется пользовательская таблица признаков: {custom_table}")
        return custom_table
    
    default_path = os.path.join("tables", "feature_table.tsv")
    
    if os.path.exists(default_path):
        if show_progress:
            print(f"Найдена существующая таблица признаков: {default_path}")
        return default_path
    
    return download_feature_table(reference_id, show_progress=show_progress)


def merge_vcfs(vcf_files, merged_output):
    """
    Объединяет VCF файлы, группируя записи по (CHROM, POS, REF),
    объединяет ALT (без дублирования) и выбирает максимальное значение QUAL.
    """
    # Словарь для группировки записей: ключ – (CHROM, POS, REF)
    records = defaultdict(list)
    # Для записи выходного VCF используем заголовок из первого входного файла.
    first_vcf = VCF(vcf_files[0])
    writer = Writer(merged_output, first_vcf)

    # Обрабатываем каждый входной VCF файл
    for vcf_file in vcf_files:
        vcf = VCF(vcf_file)
        for rec in vcf:
            key = (rec.CHROM, rec.POS, rec.REF)
            records[key].append(rec)

    # Объединяем записи по ключу
    for key, rec_list in records.items():
        if len(rec_list) == 1:
            writer.write_record(rec_list[0])
        else:
            # Объединяем ALT, убирая дублирование, и выбираем максимальное значение QUAL
            alt_set = set()
            qual_list = []
            for rec in rec_list:
                for alt in rec.ALT:
                    alt_set.add(alt)
                qual_list.append(rec.QUAL)
            merged = rec_list[0]
            merged.ALT = list(alt_set)
            merged.QUAL = max(qual_list)
            # Дополнительная логика для объединения INFO может быть добавлена здесь
            writer.write_record(merged)

    writer.close()


def run_command(command):
    return subprocess.run(command, check=True)


def check_files(vcf_files):
    """Проверяет наличие сжатых файлов и их индексов (.vcf.gz и соответствующего .tbi)"""
    for vcf in vcf_files:
        if not vcf.endswith(".vcf.gz"):
            raise ValueError(f"Файл {vcf} должен иметь расширение .vcf.gz")
        if not os.path.exists(vcf):
            raise FileNotFoundError(f"Файл {vcf} не найден")
        tbi_file = vcf + ".tbi"
        if not os.path.exists(tbi_file):
            raise FileNotFoundError(f"Индекс {tbi_file} не найден")

def get_chromosome_names_from_database(reference_id, snpeff_data_path):
    """
    Определяет имена хромосом из файла sequences.fa в базе данных.
    
    Возвращает:
    - словарь соответствия всех хромосом: {original: database_name}
    - имя первой хромосомы (для случая, когда нужно одно имя)
    """
    sequences_file = os.path.join(snpeff_data_path, reference_id, "sequences.fa")
    chromosome_names = []
    
    if os.path.exists(sequences_file):
        try:
            with open(sequences_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        # Формат FASTA: >chromosome_name [optional description]
                        # Извлекаем только имя хромосомы (до первого пробела)
                        chrom_name = line[1:].split()[0]
                        chromosome_names.append(chrom_name)
            
            if chromosome_names:
                print(f"Найдены следующие хромосомы в базе данных: {', '.join(chromosome_names)}")
                return chromosome_names, chromosome_names[0]
        except Exception as e:
            print(f"Ошибка при чтении файла последовательностей: {str(e)}")
    
    print("Не удалось определить имена хромосом из базы данных. Используется 'Chromosome' по умолчанию.")
    return [], "Chromosome"  # Значения по умолчанию

def rename_chromosomes(input_vcf, output_vcf, new_chrs="Chromosome", chrom_mapping=None):
    """
    Переименование хромосом в VCF файле.
    
    Параметры:
    - input_vcf: входной VCF файл
    - output_vcf: выходной VCF файл
    - new_chrs: имя хромосомы по умолчанию, если не задано сопоставление
    - chrom_mapping: словарь соответствия {vcf_chrom_name: db_chrom_name}
    """
    # Создание временного файла для переименования
    with tempfile.NamedTemporaryFile(mode='w+t', delete=False) as tmp_file:
        tmp_filename = tmp_file.name
        
        # Извлечение всех уникальных хромосом и замена на новое имя
        try:
            chrom_output = subprocess.check_output(
                ["bcftools", "query", "-f", '%CHROM\n', input_vcf],
                universal_newlines=True
            )
            chroms = set(chrom_output.strip().split('\n'))
            
            # Запись в формате <old_name><tab><new_name>
            for chrom in chroms:
                if not chrom:  # Пропускаем пустые строки
                    continue
                    
                # Определяем новое имя хромосомы
                if chrom_mapping and chrom in chrom_mapping:
                    # Если есть сопоставление для этой хромосомы
                    target_name = chrom_mapping[chrom]
                    print(f"Переименование хромосомы: {chrom} -> {target_name}")
                else:
                    # Используем имя по умолчанию
                    target_name = new_chrs
                    print(f"Переименование хромосомы: {chrom} -> {target_name} (по умолчанию)")
                
                tmp_file.write(f"{chrom}\t{target_name}\n")
            
            tmp_file.flush()
        except subprocess.CalledProcessError as e:
            print(f"Ошибка при получении хромосом: {e}")
            os.unlink(tmp_filename)
            raise

    # Выполнение переименования
    try:
        subprocess.run(
            ["bcftools", "annotate", "--rename-chrs", tmp_filename, "-o", output_vcf, "-O", "z", input_vcf],
            check=True
        )
    finally:
        # Удаление временного файла
        os.unlink(tmp_filename)

def annotate_vcf(input_vcf, output_vcf, reference_id):
    """
    Аннотирование VCF файла с помощью snpEff.
    
    Приоритеты:
    1. Использовать локальную базу данных, если она существует
    2. Если локальной базы нет, позволить snpEff попытаться скачать её
    3. Если база недоступна, выдать информативную ошибку
    """
    try:
        # Поиск директории с данными snpEff в conda-окружении
        conda_prefix = os.environ.get('CONDA_PREFIX')
        snpeff_data_dir = None
        
        if conda_prefix:
            # Ищем директорию данных snpEff в conda-окружении
            possible_paths = [
                os.path.join(conda_prefix, 'share', 'snpeff', 'data'),
                os.path.join(conda_prefix, 'share', 'snpeff-*', 'data'),
                os.path.join(conda_prefix, 'data')
            ]
            
            for path_pattern in possible_paths:
                matching_paths = glob.glob(path_pattern)
                for path in matching_paths:
                    if os.path.exists(path) and os.path.isdir(path):
                        snpeff_data_dir = path
                        break
                if snpeff_data_dir:
                    break
            
            if not snpeff_data_dir:
                # Если не нашли в стандартных местах, попробуем найти через запуск which
                result = subprocess.run(["which", "snpEff"], capture_output=True, text=True)
                if result.returncode == 0:
                    snpeff_path = result.stdout.strip()
                    if os.path.exists(snpeff_path):
                        with open(snpeff_path, 'r') as f:
                            content = f.read()
                            for line in content.splitlines():
                                if 'data_dir=' in line:
                                    match = re.search(r'data_dir="([^"]+)"', line)
                                    if match:
                                        snpeff_data_dir = match.group(1)
                                        break
        
        # Если всё ещё не нашли, проверяем стандартные места
        if not snpeff_data_dir:
            home = os.path.expanduser("~")
            common_paths = [
                os.path.join(home, 'snpEff', 'data'),
                '/usr/share/snpEff/data',
                '/usr/local/share/snpEff/data'
            ]
            for path in common_paths:
                if os.path.exists(path) and os.path.isdir(path):
                    snpeff_data_dir = path
                    break
        
        # Базовая команда snpEff
        cmd = ["snpEff", "ann", "-noLog", "-noStats", "-no-downstream", "-no-upstream",
               "-no-utr", "-o", "vcf"]
        
        # Проверка наличия локальной базы данных
        local_db_exists = False
        
        if snpeff_data_dir:
            print(f"Использую директорию данных snpEff: {snpeff_data_dir}")
            snpeff_home = os.path.dirname(snpeff_data_dir)
            snpeff_data_path = os.path.join(snpeff_home, "data")
            
            # Проверка существования директории базы данных
            ref_db_path = os.path.join(snpeff_data_path, reference_id)
            if os.path.exists(ref_db_path) and os.path.isdir(ref_db_path):
                print(f"Найдена локальная база данных для {reference_id} в {ref_db_path}")
                local_db_exists = True
                
                # Проверка наличия основного файла предиктора
                predictor_file = os.path.join(ref_db_path, "snpEffectPredictor.bin")
                if not os.path.exists(predictor_file):
                    print(f"Внимание: База данных существует, но отсутствует файл snpEffectPredictor.bin")
                    print(f"Возможно, база данных не была полностью построена. Рекомендуется выполнить:")
                    print(f"cd {snpeff_home}")
                    print(f"snpEff build -dataDir data -gff3 -v {reference_id}")
                    local_db_exists = False
            else:
                print(f"Локальная база данных для {reference_id} не найдена в {ref_db_path}")
                
                # Выводим доступные базы данных для информации
                available_dbs = [d for d in os.listdir(snpeff_data_path) 
                               if os.path.isdir(os.path.join(snpeff_data_path, d))]
                if available_dbs:
                    print(f"Доступные локальные базы данных:")
                    for db in available_dbs[:10]:
                        print(f"  - {db}")
                    if len(available_dbs) > 10:
                        print(f"  ... и еще {len(available_dbs) - 10} баз данных")
            
            # Добавляем путь к данным в команду
            cmd.extend(["-dataDir", snpeff_data_path])
        else:
            print("Внимание: Не удалось найти директорию данных snpEff")
        
        # Если база данных существует локально, отключаем загрузку из интернета
        if local_db_exists:
            cmd.append("-nodownload")
            print("Используем локальную базу данных без загрузки из интернета")
        else:
            print("Локальная база данных не найдена. snpEff попытается скачать базу данных...")
            
        # Добавляем имя референса и входной файл
        cmd.extend([reference_id, input_vcf])
        
        # Запускаем команду аннотации
        print(f"Запуск команды: {' '.join(cmd)}")
        with open(output_vcf, "w") as output_file:
            result = subprocess.run(cmd, stdout=output_file, stderr=subprocess.PIPE, text=True, check=False)
            
            # Проверяем результат выполнения
            if result.returncode != 0:
                stderr_output = result.stderr
                print(f"Ошибка при выполнении snpEff:")
                print(stderr_output)
                
                # Проверяем наличие характерных сообщений об ошибках
                if "Cannot find" in stderr_output and reference_id in stderr_output:
                    print(f"\nОШИБКА: База данных '{reference_id}' не найдена локально и не может быть загружена.")
                    print(f"Для использования этой базы данных необходимо сначала создать её:")
                    if snpeff_home:
                        print(f"1. Создайте директории:")
                        print(f"   mkdir -p {snpeff_data_path}/{reference_id}")
                        print(f"2. Скопируйте файлы генома и аннотаций:")
                        print(f"   cp genome.fa {snpeff_data_path}/{reference_id}/sequences.fa")
                        print(f"   cp annotations.gff {snpeff_data_path}/{reference_id}/genes.gff")
                        print(f"3. Добавьте запись в конфигурационный файл:")
                        print(f"   echo '{reference_id}.genome : {reference_id}' >> {snpeff_home}/snpEff.config")
                        print(f"4. Постройте базу данных:")
                        print(f"   cd {snpeff_home}")
                        print(f"   snpEff build -dataDir data -gff3 -v {reference_id}")
                
                # Прокидываем ошибку дальше
                result.check_returncode()  # Это вызовет исключение с кодом возврата
    
    except Exception as e:
        print(f"Ошибка при аннотации: {str(e)}")
        raise

def prepare_vcf(vcf_files, output_vcf, reference_id, show_progress=False):
    """
    Объединяет подготовку VCF-файлов:
    1. Объединяет VCF файлы
    2. Переименовывает хромосомы, автоматически определяя имена из базы данных
    3. Аннотирует VCF
    """
    try:
        # Проверка входных файлов
        check_files(vcf_files)
        
        # Создаем временные файлы для промежуточных результатов
        with tempfile.NamedTemporaryFile(suffix='.vcf.gz', delete=False) as merged_file, \
             tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as renamed_file:
            
            merged_filename = merged_file.name
            renamed_filename = renamed_file.name
        
        if show_progress:
            print("Объединяем VCF файлы...")
        merge_vcfs(vcf_files, merged_filename)
        
        # Находим директорию данных snpEff
        conda_prefix = os.environ.get('CONDA_PREFIX')
        snpeff_data_dir = None
        
        # Ищем директорию данных (этот код дублирует часть функции annotate_vcf)
        if conda_prefix:
            possible_paths = [
                os.path.join(conda_prefix, 'share', 'snpeff', 'data'),
                os.path.join(conda_prefix, 'share', 'snpeff-*', 'data'),
                os.path.join(conda_prefix, 'data')
            ]
            
            for path_pattern in possible_paths:
                matching_paths = glob.glob(path_pattern)
                for path in matching_paths:
                    if os.path.exists(path) and os.path.isdir(path):
                        snpeff_data_dir = path
                        break
                if snpeff_data_dir:
                    break
        
        # Если все еще не нашли, ищем в стандартных местах
        if not snpeff_data_dir:
            home = os.path.expanduser("~")
            common_paths = [
                os.path.join(home, 'snpEff', 'data'),
                '/usr/share/snpEff/data',
                '/usr/local/share/snpEff/data'
            ]
            for path in common_paths:
                if os.path.exists(path) and os.path.isdir(path):
                    snpeff_data_dir = path
                    break
        
        # Определяем имена хромосом из базы данных
        if snpeff_data_dir:
            snpeff_home = os.path.dirname(snpeff_data_dir)
            snpeff_data_path = os.path.join(snpeff_home, "data")
            
            if show_progress:
                print(f"Определяем имена хромосом из базы данных {reference_id}...")
            
            # Получаем список хромосом из базы данных
            chromosome_names, default_chrom = get_chromosome_names_from_database(
                reference_id, snpeff_data_path
            )
            
            # Теперь получаем хромосомы из VCF файла для сопоставления
            try:
                vcf_chroms = set()
                chrom_output = subprocess.check_output(
                    ["bcftools", "query", "-f", '%CHROM\n', merged_filename],
                    universal_newlines=True
                )
                for chrom in chrom_output.strip().split('\n'):
                    if chrom:
                        vcf_chroms.add(chrom)
                
                if show_progress:
                    print(f"Хромосомы в VCF файле: {', '.join(vcf_chroms)}")
                
                # Создаем сопоставление имен хромосом
                # Если число хромосом совпадает, делаем прямое сопоставление
                chrom_mapping = {}
                
                if len(vcf_chroms) == len(chromosome_names):
                    # Простое сопоставление по порядку
                    for vcf_chrom, db_chrom in zip(sorted(vcf_chroms), sorted(chromosome_names)):
                        chrom_mapping[vcf_chrom] = db_chrom
                    
                    if show_progress:
                        print("Созданы сопоставления хромосом на основе порядкового номера:")
                        for vcf_chrom, db_chrom in chrom_mapping.items():
                            print(f"  {vcf_chrom} -> {db_chrom}")
                else:
                    # Если число хромосом не совпадает, используем имя первой хромосомы для всех
                    if chromosome_names:
                        for vcf_chrom in vcf_chroms:
                            chrom_mapping[vcf_chrom] = default_chrom
                        
                        if show_progress:
                            print(f"Все хромосомы будут переименованы в {default_chrom}")
            
            except Exception as e:
                print(f"Ошибка при определении хромосом из VCF: {str(e)}")
                chrom_mapping = None
        else:
            chrom_mapping = None
            default_chrom = "Chromosome"
        
        # Переименование хромосом с использованием полученного сопоставления
        if show_progress:
            print("Переименовываем хромосомы...")
        
        rename_chromosomes(merged_filename, renamed_filename, default_chrom, chrom_mapping)
        
        if show_progress:
            print(f"Аннотируем VCF файл с использованием референса {reference_id}...")
        annotate_vcf(renamed_filename, output_vcf, reference_id)
        
        # Удаляем временные файлы
        os.unlink(merged_filename)
        os.unlink(renamed_filename)
        
        if show_progress:
            print(f"VCF файл подготовлен: {output_vcf}")
        
        return output_vcf
        
    except Exception as e:
        print(f"Ошибка при подготовке VCF: {str(e)}")
        sys.exit(1)


def save_with_progress(df, filename, show_progress=True, chunk_size=1000):
    """Оптимизированное сохранение с прогресс-баром"""
    # Создаем книгу и лист
    wb = Workbook()
    ws = wb.active
    
    # Записываем заголовки
    ws.append(list(df.columns))
    
    # Конвертируем DataFrame в список списков для быстрого доступа
    data = df.values.tolist()
    
    # Настройка прогресс-бара
    if show_progress:
        pbar = tqdm(total=len(data), desc="Сохранение в Excel", unit="row")
    
    # Запись чанками
    for i in range(0, len(data), chunk_size):
        chunk = data[i:i+chunk_size]
        
        # Добавляем чанк данных
        for row in chunk:
            ws.append(row)
        
        # Обновляем прогресс
        if show_progress:
            pbar.update(len(chunk))
    
    # Финализация
    if show_progress:
        pbar.close()
    
    wb.save(filename)


def read_vcf_with_pysam(file_path):
    """
    Чтение VCF файла с использованием pysam и преобразование в DataFrame с колонками POS, REF, Allele.
    """
    data = []
    with pysam.VariantFile(file_path) as vcf:
        for record in vcf.fetch():
            pos = record.pos
            ref = record.ref
            alts = ",".join(record.alts) if record.alts else ""
            data.append([pos, ref, alts])

    return pd.DataFrame(data, columns=["POS", "REF", "Allele"])


def extract_alt(alt_value):
    """Извлечение значения из строки формата 'value=...'."""
    match = re.search(r"value='(.*?)'", alt_value)
    return match.group(1) if match else None


def load_feature_tables(file_path):
    """Загрузка таблиц признаков."""
    df = pd.read_csv(file_path, sep='\t')
    feature_table = df[['symbol', 'locus_tag']].drop_duplicates().dropna()
    name_table = df[['name', 'locus_tag']].drop_duplicates()
    return feature_table, name_table


def get_locus_tag(symbol_value, feature_table):
    """Получение locus_tag по символу."""
    if not symbol_value or symbol_value == "":
        return "None"
        
    locus_tags = feature_table[feature_table['symbol'] == symbol_value]['locus_tag'].values
    return locus_tags[0] if locus_tags.size > 0 else symbol_value


def get_name_tag(locus_tag_value, name_table):
    """Получение name по locus_tag с улучшенной обработкой ошибок."""
    if not locus_tag_value or locus_tag_value == "":
        return "None"
        
    locus_names = locus_tag_value.split("-")
    name_tag = []
    
    for locus_name in locus_names:
        if not locus_name:
            name_tag.append("None")
            continue
            
        name_entries = name_table[name_table['locus_tag'] == locus_name]['name'].values
        
        if len(name_entries) == 0:
            name_tag.append("None")
        elif len(name_entries) == 1:
            name = str(name_entries[0]) if name_entries[0] is not None else "None"
            name_tag.append(name)
        else:
            name = str(name_entries[1]) if len(name_entries) > 1 and name_entries[1] is not None else "None"
            name_tag.append(name)
    
    return ",".join(name_tag)


def rename_gene_id(gene_name, gene_id, feature_table):
    """Обновление gene_id на основе feature_table с улучшенной обработкой ошибок."""
    # Проверка на пустые значения
    if not gene_name or not gene_id:
        return gene_id if gene_id else "None"
    
    gene_id_parts = gene_id.split("-")
    gene_name_parts = gene_name.split("-")
    
    # Убедимся, что gene_name_parts имеет как минимум столько же элементов, сколько gene_id_parts
    while len(gene_name_parts) < len(gene_id_parts):
        gene_name_parts.append("")
    
    # Обновляем каждую часть gene_name на соответствующий locus_tag
    for index, part in enumerate(gene_id_parts):
        if index < len(gene_name_parts):  # Безопасная проверка
            gene_name_parts[index] = get_locus_tag(gene_name_parts[index], feature_table)
    
    return "-".join(gene_name_parts)


def process_vcf_record(record, feature_table, name_table):
    """Обработка одной записи VCF с улучшенной обработкой ошибок."""
    try:
        # Безопасно получаем annotations
        annotations = record.info.get('ANN', [])
        rows = []
        
        # Инициализируем annotation_parts с пустыми значениями
        annotation_parts = [""]*16
        
        if annotations and len(annotations) > 0:
            first_annotation = annotations[0].split("|")
            # Безопасно заполняем annotation_parts
            for i in range(min(len(first_annotation), 16)):
                annotation_parts[i] = first_annotation[i]
        
        # Безопасно получаем allele
        allele = ",".join(record.alts) if record.alts else ""
        
        # Обработка аннотаций с проверкой на наличие alts
        if record.alts and len(record.alts) == 1:
            # Безопасно получаем данные для одного альтернативного аллеля
            temp_hgvsc = annotation_parts[9] if len(annotation_parts) > 9 else ""
            temp_hgvsp = annotation_parts[10] if len(annotation_parts) > 10 else ""
            temp_annotation = annotation_parts[1] if len(annotation_parts) > 1 else ""
            temp_putative_impact = annotation_parts[2] if len(annotation_parts) > 2 else ""
            temp_feature_type = annotation_parts[5] if len(annotation_parts) > 5 else ""
        else:
            # Инициализируем списки для нескольких аннотаций
            temp_hgvsc = []
            temp_hgvsp = []
            temp_annotation = []
            temp_putative_impact = []
            temp_feature_type = []
            
            for ann in annotations:
                ann_split = ann.split("|")
                # Безопасно добавляем данные в списки
                temp_annotation.append(ann_split[1] if len(ann_split) > 1 else "")
                temp_putative_impact.append(ann_split[2] if len(ann_split) > 2 else "")
                temp_feature_type.append(ann_split[5] if len(ann_split) > 5 else "")
                temp_hgvsc.append(ann_split[9] if len(ann_split) > 9 else "")
                temp_hgvsp.append(ann_split[10] if len(ann_split) > 10 else "")
            
            # Объединяем списки в строки
            temp_hgvsp = ",".join(temp_hgvsp)
            temp_hgvsc = ",".join(temp_hgvsc)
            temp_annotation = ",".join(temp_annotation)
            temp_putative_impact = ",".join(temp_putative_impact)
            temp_feature_type = ",".join(temp_feature_type)
        
        # Безопасно получаем gene_name и gene_id
        gene_name = annotation_parts[3] if len(annotation_parts) > 3 else ""
        gene_id_raw = annotation_parts[4] if len(annotation_parts) > 4 else ""
        
        # Безопасное переименование gene_id
        try:
            gene_id = rename_gene_id(gene_name, gene_id_raw, feature_table)
        except Exception as e:
            print(f"Ошибка при переименовании gene_id: {e}")
            gene_id = gene_id_raw
        
        # Безопасное получение name_tag
        try:
            name_tag = get_name_tag(gene_id, name_table)
        except Exception as e:
            print(f"Ошибка при получении name_tag: {e}")
            name_tag = "None"
        
        # Формируем строку данных с безопасным доступом к элементам
        row = [
            record.pos,
            record.ref,
            allele,
            temp_annotation,
            temp_putative_impact,
            gene_name,
            gene_id,
            name_tag,
            temp_feature_type,
            annotation_parts[7] if len(annotation_parts) > 7 else "",
            temp_hgvsc,
            temp_hgvsp,
            annotation_parts[11] if len(annotation_parts) > 11 else "",
            annotation_parts[12] if len(annotation_parts) > 12 else "",
            annotation_parts[13] if len(annotation_parts) > 13 else "",
            annotation_parts[15] if len(annotation_parts) > 15 else ""
        ]
        
        rows.append(row)
        return rows
    
    except Exception as e:
        print(f"Ошибка при обработке записи VCF на позиции {record.pos}: {e}")
        # Возвращаем запись с минимальной информацией в случае ошибки
        return [[record.pos, record.ref, "", "", "", "", "", "None", "", "", "", "", "", "", "", "Error"]]


def process_vcf_parallel(input_vcf, feature_table, name_table, show_progress, threads=4):
    """Параллельная обработка VCF файла с сохранением порядка записей."""
    columns = [
        "POS", "REF", "Allele", "Annotation", "Putative_impact", "Gene Name", "Gene ID", "name",
        "Feature type", "Transcript biotype", "HGVS.c", "HGVS.p", "cDNA_position / cDNA_len",
        "CDS_position / CDS_len", "Protein_position / Protein_len", "Errors"
    ]

    with pysam.VariantFile(input_vcf) as vcf:
        records = list(vcf.fetch())

    process_record = partial(process_vcf_record, 
                           feature_table=feature_table, 
                           name_table=name_table)

    data = []
    with ThreadPoolExecutor(max_workers=threads) as executor:
        if show_progress:
            records_iter = tqdm(records, desc="Обработка записей", unit="rec")
        else:
            records_iter = records

        results = executor.map(process_record, records_iter)
        
        for result in results:
            data.extend(result)

    return pd.DataFrame(data, columns=columns).sort_values(by="POS")


def extract_annotations(input_vcf, output_excel, show_progress, threads, feature_table_path):
    """Основная функция для извлечения аннотаций из VCF."""
    try:
        feature_table, name_table = load_feature_tables(feature_table_path)
        result_df = process_vcf_parallel(
            input_vcf, feature_table, name_table, show_progress=show_progress, threads=threads
        )
        
        if len(result_df) > 100_000:
            result_df.to_excel(output_excel, index=False, engine='openpyxl')
        else:
            save_with_progress(result_df, output_excel, show_progress)
            
        if show_progress:
            print(f"Файл {output_excel} был успешно создан!")
        
        return result_df
    except Exception as e:
        print(f"Ошибка при извлечении аннотаций: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def update_excel_with_vcfs(vcfs, input_file, output_file, show_progress, threads):
    """Обновление Excel файла данными из VCF с сохранением порядка файлов."""
    try:
        df = pd.read_excel(input_file).set_index('POS')
    
        def safe_process_vcf(vcf):
            try:
                df_vcf = read_vcf_with_pysam(vcf)
                file_name = os.path.basename(vcf).split(".")[0]
                return (file_name, df_vcf.set_index('POS')['Allele'])
            except Exception as e:
                print(f"Ошибка обработки {vcf}: {e}")
                return None
    
        with ThreadPoolExecutor(max_workers=threads) as executor:
            if show_progress:
                vcfs = tqdm(vcfs, desc="Обработка VCF", unit="file")
    
            results = list(executor.map(safe_process_vcf, vcfs))
    
        for result in filter(None, results):
            file_name, series = result
            df[file_name] = series
    
        df.reset_index(inplace=True)
    
        df.to_excel(output_file, index=False, engine='openpyxl')
        
        if show_progress:
            print(f"Файл сохранен как {output_file}")
        
        return df
    except Exception as e:
        print(f"Ошибка при обновлении Excel: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def process_all(vcf_files, output_excel, reference_id, custom_table=None, temp_vcf=None, 
                show_progress=True, threads=4):
    """
    Полный процесс обработки:
    1. Подготовка и аннотирование VCF
    2. Извлечение аннотаций
    3. Обновление с данными из всех VCF
    """
    # Если временный VCF не указан, создаем его
    if temp_vcf is None:
        temp_vcf = "annotated_merged.vcf"
    
    # Создаем временный Excel
    temp_excel = "temp_extracted.xlsx"
    
    try:
        # Обеспечиваем наличие таблицы признаков
        feature_table_path = ensure_feature_table(reference_id, custom_table, show_progress)
        
        # Шаг 1: Подготовка VCF
        if show_progress:
            print("=== Шаг 1: Подготовка и аннотирование VCF ===")
        prepare_vcf(vcf_files, temp_vcf, reference_id, show_progress)
        
        # Шаг 2: Извлечение аннотаций
        if show_progress:
            print("\n=== Шаг 2: Извлечение аннотаций ===")
        extract_annotations(temp_vcf, temp_excel, show_progress, threads, feature_table_path)
        
        # Шаг 3: Обновление с данными из всех VCF
        if show_progress:
            print("\n=== Шаг 3: Обновление с данными из всех VCF ===")
        update_excel_with_vcfs(vcf_files, temp_excel, output_excel, show_progress, threads)
        
        # Удаление временных файлов
        if os.path.exists(temp_vcf):
            os.remove(temp_vcf)
        if os.path.exists(temp_excel):
            os.remove(temp_excel)
        
        if show_progress:
            print(f"\nПроцесс успешно завершен! Результат: {output_excel}")
            
    except Exception as e:
        print(f"Ошибка в процессе обработки: {str(e)}")
        import traceback
        traceback.print_exc()
        # Удаляем временные файлы в случае ошибки
        if os.path.exists(temp_vcf):
            os.remove(temp_vcf)
        if os.path.exists(temp_excel):
            os.remove(temp_excel)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Универсальный скрипт для комплексной работы с VCF и Excel.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Общие аргументы для всех команд, связанных с референсом
    ref_args = argparse.ArgumentParser(add_help=False)
    ref_args.add_argument("-r", "--reference", default="Mycobacterium_tuberculosis_h37rv",
                          help="Референсный геном для использования со snpEff (по умолчанию: Mycobacterium_tuberculosis_h37rv)")
    ref_args.add_argument("--feature-table", help="Путь к пользовательской таблице признаков TSV")
    
    # Команда extract
    extract_parser = subparsers.add_parser("extract", parents=[ref_args], help="Извлечь аннотации из VCF в Excel.")
    extract_parser.add_argument("-i", "--input", required=True, help="Путь к входному VCF файлу.")
    extract_parser.add_argument("-o", "--output", required=True, help="Путь к выходному Excel файлу.")
    extract_parser.add_argument("-p", "--progress", action="store_true", help="Отображать прогресс выполнения.")
    extract_parser.add_argument("-th", "--threads", type=int, default=4, help="Количество потоков для обработки.")

    # Команда update
    update_parser = subparsers.add_parser("update", help="Обновить Excel с помощью данных из VCF.")
    update_parser.add_argument("-v", "--vcfs", nargs='+', required=True, help="Список VCF файлов.")
    update_parser.add_argument("-t", "--table", required=True, help="Путь к Excel таблице.")
    update_parser.add_argument("-o", "--output", required=True, help="Путь к выходному Excel файлу.")
    update_parser.add_argument("-p", "--progress", action="store_true", help="Отображать прогресс выполнения.")
    update_parser.add_argument("-th", "--threads", type=int, default=4, help="Количество потоков для обработки.")

    # Команда prepare
    prepare_parser = subparsers.add_parser("prepare", parents=[ref_args], help="Подготовить и аннотировать VCF файлы.")
    prepare_parser.add_argument("-i", "--input", nargs='+', required=True, help="Список VCF файлов.")
    prepare_parser.add_argument("-o", "--output", required=True, help="Путь к выходному VCF файлу.")
    prepare_parser.add_argument("-p", "--progress", action="store_true", help="Отображать прогресс выполнения.")

    # НОВАЯ команда - all-in-one
    all_parser = subparsers.add_parser("all", parents=[ref_args],
                                       help="Выполнить все шаги: подготовить, извлечь и обновить в одной команде.")
    all_parser.add_argument("-v", "--vcfs", nargs='+', required=True, 
                           help="Список VCF файлов (.vcf.gz с индексами .tbi).")
    all_parser.add_argument("-o", "--output", required=True, 
                           help="Путь к итоговому Excel файлу.")
    all_parser.add_argument("-p", "--progress", action="store_true", 
                           help="Отображать прогресс выполнения.")
    all_parser.add_argument("-th", "--threads", type=int, default=4, 
                           help="Количество потоков для обработки.")
    all_parser.add_argument("-t", "--temp", 
                           help="Путь к временному VCF файлу (не будет удален после выполнения).")

    # Команда download-table - самостоятельная функция для скачивания таблицы признаков
    dl_parser = subparsers.add_parser("download-table", parents=[ref_args],
                                      help="Скачать таблицу признаков для указанного референса.")
    dl_parser.add_argument("-o", "--output", 
                          help="Директория для сохранения таблицы признаков (по умолчанию: tables)")
    dl_parser.add_argument("-p", "--progress", action="store_true", help="Отображать прогресс выполнения.")

    args = parser.parse_args()

    if args.command == "extract":
        feature_table_path = ensure_feature_table(args.reference, args.feature_table, args.progress)
        extract_annotations(args.input, args.output, args.progress, args.threads, feature_table_path)
    elif args.command == "update":
        update_excel_with_vcfs(args.vcfs, args.table, args.output, args.progress, args.threads)
    elif args.command == "prepare":
        prepare_vcf(args.input, args.output, args.reference, args.progress)
    elif args.command == "all":
        process_all(args.vcfs, args.output, args.reference, args.feature_table, 
                   args.temp, args.progress, args.threads)
    elif args.command == "download-table":
        output_dir = args.output if args.output else "tables"
        download_feature_table(args.reference, output_dir, args.progress)

if __name__ == "__main__":
    main()
