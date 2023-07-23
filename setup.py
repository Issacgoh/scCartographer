
from setuptools import setup, find_packages

setup(
    name='MyScRNASeqPackage',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'pandas',
        'scikit-learn',
    ],
)
