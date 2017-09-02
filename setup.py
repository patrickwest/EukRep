from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='EukRep',
    version='0.6.2',
    description='Classification of Eukaryotic and Prokaryotic sequences from metagenomic datasets',
    url='https://github.com/patrickwest/EukRep',
    author='Patrick West',
    author_email='patrickwest@berkeley.edu',
    license='MIT',
    package_data={'EukRep': ['models/*.pickle']},
    packages=['EukRep'],
    scripts=['bin/EukRep'],
    install_requires=[
          'numpy',
          'sklearn',
          'kpal',
          'biopython'
    ],
    zip_safe=False)
