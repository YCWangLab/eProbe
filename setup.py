
from setuptools import setup, find_packages

setup(
    name="eProbe",
    version="0.0.1",
    author="Zi-Hao Huang",
    author_email="zh384@cam.ac.uk",
    description="A package for probe design",
    long_description="A one-stop capture probe design toolkit for genetic diversity reconstructions from ancient environmental DNA",
    packages=find_packages(),
    install_requires=[ 
    'biopython>=1.83',
    'numpy>=1.24.4',
    'pandas>=2.0.3',
    'pysam>=0.22.0',
    'pybedtools>=0.10.0',
    'matplotlib==3.7.5',
    'seaborn==0.13.2',
    'scipy==1.10.1',
    'fastcluster==1.2.6'
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: CC BY-NC 4.0",
        "Operating System :: OS Independent",
    ],
    license="Creative Commons Attribution-NonCommercial 4.0",
    url="https://github.com/YCWangLab/eProbe",
    )