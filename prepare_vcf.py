#!/usr/bin/env python3

import subprocess
import argparse
import sys
import os
from collections import defaultdict
from cyvcf2 import VCF, Writer

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
    subprocess.run(command, check=True)

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

def main():
    parser = argparse.ArgumentParser(description="Обработка VCF файлов.")
    parser.add_argument("-i", "--input", type=str, nargs="+", required=True,
                        help="Список сжатых VCF файлов (.vcf.gz) для объединения")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Имя выходного VCF файла (результат аннотирования)")
    args = parser.parse_args()
    vcf_files = args.input
    output = args.output

    if not vcf_files:
        print("Ошибка: Не указаны VCF файлы")
        sys.exit(1)

    try:
        # Проверка входных файлов
        check_files(vcf_files)
        
        # Шаг 1. Объединение VCF файлов с использованием внутренней функции
        print("Объединяем VCF файлы")
        merge_vcfs(vcf_files, "merged.vcf.gz")

        # Шаг 2. Переименование хромосом
        print("Переименовываем хромосомы...")
        run_command(["scripts/rename_vcf", "merged.vcf.gz", "merged.renamed.vcf"])

        # Шаг 3. Аннотирование VCF файла с использованием snpEff
        print("Аннотируем VCF файл...")
        with open(output, "w") as output_file:
            subprocess.run(
                [
                    "java", "-jar", "/home/zerg/soft/src/snpEff/snpEff.jar",
                    "ann", "-noLog", "-noStats", "-no-downstream", "-no-upstream",
                    "-no-utr", "-o", "vcf", "Mycobacterium_tuberculosis_h37rv",
                    "merged.renamed.vcf"
                ],
                stdout=output_file,
                check=True
            )

        print("Процесс завершён!")

    except Exception as e:
        print(f"Ошибка: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
