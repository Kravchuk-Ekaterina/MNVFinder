# MNVFinder
BI 2022 Spring Student Project<br>
## Motivation
Multi-nucleotide variants (MNVs) are defined as clusters of two or more nearby variants existing on the same haplotype in an individual. When variants in an MNV are found within the same codon, the overall impact may differ from the functional consequences of the individual variants.<br>
![motivation](/images/motivation.jpg "motivation")<br>
Modern publicly available tools incorrectly annotate polymorphisms in the same codon, considering their contributions independently. It would be useful to create a tool to properly annotate MNVs. <br>
## Methods and Results
![MNVFinder](/images/MNVFinder.jpg "MNVFinder")<br>
MNVFinder is a tool for annotation of MNVs. It works with VEP annotations https://www.ensembl.org/info/docs/tools/vep/index.html.<br>
The tool takes a VEP annotation in vcf format as input data. It identifies SNPs within a siggle codon and gives the annotation of the found MNVs.<br>
You can find the validation in ./validation/test1/ and ./validation/test2/.<br>
## Requirements
The tool is Linux command line tool and requires Python >= 3.6. You must have installed pandas >= 1.4.1 library.
## Installation
Clon this repository:
```bash
git clone https://github.com/Kravchuk-Ekaterina/MNVFinder
```
## Usage
```bash
python path_to_MNVFinder/mnvfinder.py -h
```
```bash
usage: mnvfinder.py [-h] -i FILE -o OUTDIR

MNVFinder: the tool for annotation MNVs

optional arguments:
  -h, --help            show this help message and exit
  -i FILE, --input FILE
                        .tab file containing vep annotation
  -o OUTDIR, --outdir OUTDIR
                        the directory for the output file
```
The examples of usage can be found in ./validation/test1/labjournal.md and ./validation/test2/labjournal.ipynb.<br>

