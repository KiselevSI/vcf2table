#!.venv/bin/python
import argparse
import sys
import pandas as pd
import pysam
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import re
import gzip
import os
from functools import partial

from openpyxl.utils.dataframe import dataframe_to_rows
from openpyxl import Workbook

def save_with_row_progress(df, filename):
    """Запись с прогрессом по строкам (для больших файлов)"""
    with pd.ExcelWriter(filename, engine='openpyxl') as writer, \
         tqdm(total=len(df), desc="Запись строк") as pbar:
        
        df.to_excel(writer, index=False)
        sheet = writer.sheets['Sheet1']
        
        # Прогресс для форматов и стилей
        for row in range(2, len(df)+2):
            sheet.row_dimensions[row].height = 15
            pbar.update(1)
            
            
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
    locus_tags = feature_table[feature_table['symbol'] == symbol_value]['locus_tag'].values
    return locus_tags[0] if locus_tags.size > 0 else symbol_value

def get_name_tag(locus_tag_value, name_table):
    locus_names = locus_tag_value.split("-")
    name_tag = []
    for locus_name in locus_names:
        name_entries = name_table[name_table['locus_tag'] == locus_name]['name'].values
        if len(name_entries) <= 1:
            name_tag.append("None")
            continue
        name = name_entries[1] if name_entries[1] else "None"
        name_tag.append(name)
    return ",".join(name_tag)

def rename_gene_id(gene_name, gene_id, feature_table):
    """Обновление gene_id на основе feature_table."""
    gene_id_parts = gene_id.split("-")
    gene_name_parts = gene_name.split("-")

    for index, part in enumerate(gene_id_parts):
        gene_name_parts[index] = get_locus_tag(gene_name_parts[index], feature_table)

    return "-".join(gene_name_parts)

def process_vcf_record(record, feature_table, name_table):
    """Обработка одной записи VCF."""
    annotations = record.info.get('ANN', [])
    rows = []
    
    if annotations:
        first_annotation = annotations[0].split("|")
        annotation_parts = first_annotation
    else:
        annotation_parts = [""]*16


    allele = ",".join(record.alts) if record.alts else ""

    
    if len(record.alts) == 1:
        temp_hgvsc = annotation_parts[9]
        temp_hgvsp = annotation_parts[10]
        temp_annotation = annotation_parts[1]
        temp_putative_impact = annotation_parts[2]
        temp_feature_type = annotation_parts[5]

    else:
        temp_hgvsc = []
        temp_hgvsp = []
        temp_annotation = []
        temp_putative_impact = []
        temp_feature_type = []
        for ann in annotations:
            ann_split = ann.split("|")
            temp_annotation.append(ann_split[1])
            temp_putative_impact.append(ann_split[2])
            temp_feature_type.append(annotation_parts[5])
            
            temp_hgvsc.append(ann_split[9])
            temp_hgvsp.append(ann_split[10])
        temp_hgvsp = ",".join(temp_hgvsp)
        temp_hgvsc = ",".join(temp_hgvsc)
        temp_annotation = ",".join(temp_annotation)
        temp_putative_impact = ",".join(temp_putative_impact)
        temp_feature_type = ",".join(temp_feature_type)
        if len(temp_hgvsp) <= 1:
            temp_hgvsp = ""

    gene_id = rename_gene_id(
        annotation_parts[3],
        annotation_parts[4],
        feature_table
    )

    row = [
        record.pos,
        record.ref,
        allele,
        temp_annotation,
        temp_putative_impact,
        annotation_parts[3],
        gene_id,
        get_name_tag(gene_id, name_table),
        temp_feature_type,
        annotation_parts[7],
        temp_hgvsc,
        temp_hgvsp,
        annotation_parts[11],
        annotation_parts[12],
        annotation_parts[13],
        annotation_parts[15]
    ]
    
    rows.append(row)
    return rows

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
            records_iter = tqdm(records, desc="Processing records", unit="rec")
        else:
            records_iter = records

        results = executor.map(process_record, records_iter)
        
        for result in results:
            data.extend(result)

    return pd.DataFrame(data, columns=columns).sort_values(by="POS")

def extract_annotations(input_vcf, output_excel, show_progress, threads):
    """Основная функция для извлечения аннотаций из VCF."""
    annotation_table_path = "tables/feature_table.tsv"
    feature_table, name_table = load_feature_tables(annotation_table_path)
    result_df = process_vcf_parallel(
        input_vcf, feature_table, name_table, show_progress=show_progress, threads=threads
    )
    
    #if show_progress:
    #    print()  # Добавляем отступ между прогресс-барами

    #if show_progress:
    #    print("Сохранение в Excel...")
    #    save_with_row_progress(result_df, output_excel)  
    #else:  
    #    result_df.reset_index().to_excel(output_excel, index=False, engine='openpyxl')
    
    
    if len(result_df) > 100_000:
        result_df.to_excel(output_excel, index=False, engine='openpyxl')
    else:
        save_with_progress(result_df, output_excel, show_progress)
        
    
    print(f"Файл {output_excel} был успешно создан!")



def update_excel_with_vcfs(vcfs, input_file, output_file, show_progress, threads):
    """Обновление Excel файла данными из VCF с сохранением порядка файлов."""
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

    #print(df[])


    df.to_excel(output_file, index=False, engine='openpyxl')
    #if len(df) > 100_000:
    #    df.to_excel(output_file, index=False, engine='openpyxl')
    #else:
    #    save_with_progress(df, output_file, show_progress)
#
# Быстрый CSV с прогрессом
#with tqdm(total=len(df)) as pbar:
#    df.to_csv(filename, chunksize=10_000, callback=lambda _: pbar.update(10_000))
#
    print(f"Файл сохранен как {output_file}")

def main():
    parser = argparse.ArgumentParser(description="Оптимизированный скрипт для работы с VCF и Excel.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    extract_parser = subparsers.add_parser("extract", help="Извлечь аннотации из VCF в Excel.")
    extract_parser.add_argument("-i", "--input", required=True, help="Путь к входному VCF файлу.")
    extract_parser.add_argument("-o", "--output", required=True, help="Путь к выходному Excel файлу.")
    extract_parser.add_argument("-p", "--progress", action="store_true", help="Отображать прогресс выполнения.")
    extract_parser.add_argument("-th", "--threads", type=int, default=4, help="Количество потоков для обработки.")

    update_parser = subparsers.add_parser("update", help="Обновить Excel с помощью данных из VCF.")
    update_parser.add_argument("-v", "--vcfs", nargs='+', required=True, help="Список VCF файлов.")
    update_parser.add_argument("-t", "--table", required=True, help="Путь к Excel таблице.")
    update_parser.add_argument("-o", "--output", required=True, help="Путь к выходному Excel файлу.")
    update_parser.add_argument("-p", "--progress", action="store_true", help="Отображать прогресс выполнения.")
    update_parser.add_argument("-th", "--threads", type=int, default=4, help="Количество потоков для обработки.")

    args = parser.parse_args()

    if args.command == "extract":
        extract_annotations(args.input, args.output, args.progress, args.threads)
    elif args.command == "update":
        update_excel_with_vcfs(args.vcfs, args.table, args.output, args.progress, args.threads)

if __name__ == "__main__":
    main()