
from setuptools import setup, find_packages

setup(
    name="eProbe",
    version="0.0.1",
    author="Zi-Hao Huang",
    author_email="zh384@cam.ac.uk",
    description="A package for probe design",
    long_description="A one-stop target genome capture baits design toolkit for ancient environmental DNA",
    packages=find_packages(),
    install_requires=[ 
    'biopython>=1.7.0',
    'numpy>=1.19.5',
    'pandas>=1.1.5',
    'pysam>=0.11.2.2',
    'PyVCF>=0.6.8'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    license="GNU General Public License version 3",
    url="https://github.com/YCWangLab/eProbes",
    )