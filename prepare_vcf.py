#!.venv/bin/python

import subprocess
import argparse
import sys
import os

def run_command(command):
    subprocess.run(command, check=True)

def check_files(vcf_files):
    """Проверяет наличие сжатых файлов и их индексов"""
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
                        help="Имя выходного VCF файла")
    
    args = parser.parse_args()
    vcf_files = args.input
    output = args.output

    if not vcf_files:
        print("Ошибка: Не указаны VCF файлы")
        sys.exit(1)

    try:
        # Проверка входных файлов
        check_files(vcf_files)
        # bcftools merge -m all *.vcf.gz > merged.last.vcf
        # vcf-merge *vcf.gz > last.vcf
        # Объединение файлов
        print("Объединяем VCF файлы")
        merge_command = [
            'vcf-merge'
        ] + vcf_files

        with open("merged.vcf", "w") as outfile:
            subprocess.run(merge_command, stdout=outfile, check=True)

        

        # Переименование хромосом
        print("Переименовываем хромосомы...")
        run_command(["scripts/rename_vcf", "merged.vcf", "merged.renamed.vcf"])

        # Аннотирование
        print("Аннотируем VCF файл...")
        with open(output, "w") as output_file:
            subprocess.run(
                [
                    "java", "-jar", "/home/humdrum/soft/src/snpEff/snpEff.jar",
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