# -*- coding: utf-8 -*-
import io
import re
from setuptools import setup, find_packages
import sys

# BU SATIRLAR SORUNUN KALICI ÇÖZÜMÜDÜR.
# Python'a, README.md dosyasını hangi işletim sisteminde olursa olsun
# her zaman UTF-8 kodlamasıyla okumasını söylüyoruz.
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

def get_version():
    with open('mlfp/__init__.py', 'r', encoding='utf-8') as f:
        content = f.read()
    match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", content, re.M)
    if match:
        return match.group(1)
    raise RuntimeError("Unable to find version string.")

setup(
    name="mlfp",
    version=get_version(),
    description="mlfp: Mittag-Leffler and Prabhakar function, MLF",
    long_description=long_description,
    long_description_content_type="text/markdown", # Bu satır da önemlidir
    author="Mehmet Keçeci",
    maintainer="Mehmet Keçeci",
    author_email="enfo@tuta.io",
    maintainer_email="enfo@tuta.io",
    url="https://github.com/WhiteSymmetry/mlfp",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "matplotlib",
    ],
    extras_require={
        "all": [],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT",
        "Operating System :: OS Independent"
    ],
    python_requires='>=3.11',
    license="MIT",
)

