from setuptools import setup, find_packages

def parse_requirements(filename):
    """Читает зависимости из файла requirements.txt"""
    with open(filename) as f:
        return f.read().splitlines()

setup(
    name="vcf2table",  # Уникальное имя вашего проекта
    version="3.0.0",  # Версия проекта
    description="Скрипт для обработки VCF файлов и обновления таблиц",  # Краткое описание
    long_description=open("README.md").read(),  # Длинное описание, можно использовать README.md
    long_description_content_type="text/markdown",  # Указываем, что README в формате Markdown
    author="Ваше Имя",  # Ваше имя
    author_email="your_email@example.com",  # Ваш email
    url="https://github.com/your_username/your_repository",  # Ссылка на ваш проект
    packages=find_packages(),  # Автоматически найдет все пакеты и модули
    install_requires=parse_requirements('requirements.txt'),  # Загрузка зависимостей из requirements.txt
    entry_points={
        "console_scripts": [
            'vcf2table=vcf2table:main',  # Свяжет скрипт с командой в терминале
        ]
    },
    classifiers=[
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.12",  # Указываем минимальную версию Python
)
