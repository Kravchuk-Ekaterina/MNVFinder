# Running the complete pipeline for MNV annotation
## 1. Converting the raw 23andMe data to vcf
I use plink - tool widely used in population genetics https://www.cog-genomics.org/plink/
I remove all SNPs corresponding to deletions and insertions, to make the file compatible with annotation tools
```bash
path_to_plink/plink --23file <23_and_me_input.txt> --recode vcf --out snps_clean --output-chr MT --snps-only just-acgt
```
## 2. Running VEP
```bash
path_to_vep/vep --cache -canonical -i ./data/test1/snps_clean.vcf -o ./data/test1/test1_vep_coding.vcf --coding_only
```
## 3. Running MNVFinder
```bash
python ../../code/mnvfinder.py -i test1_vep_coding.vep -o .
```
The output:
```bash
Reading data...
Searching for polymorphisms in the same codon...
Would like to save file with MNVs annotated as SNPs? y/n y
Would like to save file with SNP annotation without misannotated MNVs? y/n y
Annotating MNVs...
Succesfully finished
Found 60 MNVs:
missense_variant 59
stop_lost 1
The output file(s) can be found in  .
```
