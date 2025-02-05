
если нужно создать индексные файлы .vcf.gz.tbi
ls *.vcf.gz | parallel tabix -p vcf

Сначало надо слить все vcf файлы в 1 
bcftools merge *.vcf -o merged.vcf

изменить имя хромосов в этом vcf файле, так как snpEff выдает ошибку из-за неправильного имени хромосомы, есть флаг -h для справки
rename_vcf merged.vcf output.vcf

аннотация vcf 
snpeff ann -noLog -noStats -no-downstream -no-upstream -no-utr -o vcf Mycobacterium_tuberculosis_h37rv output.vcf > output.annotated.vcf

Создание виртуального окружения
python3 -m venv .venv

Активация виртуального окружения
source .venv/bin/activate

Установка зависимостей через setup.py
pip install -e .

Сделать исполняемым
chmod +x get_data_from_vcf.py get_table.py

Генерация главной таблицы, есть флаг -h для справки
./get_data_from_vcf.py -i merge.vcf -o output.xlsx

Добавление в главную таблицу новых вариантов, есть флаг -h для справки
./get_table.py -v vcfs/*vcf.gz -t merged.xlsx