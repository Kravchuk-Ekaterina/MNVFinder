# Importing modules

import pandas as pd
import os
import warnings
import argparse

# Parser for command line options

from argparse import ArgumentParser

parser = ArgumentParser(description="MNVFinder: the tool for annotation MNVs")
parser.add_argument("-i", "--input", dest="filename", required=True, help=".tab file containing vep annotation", metavar="FILE")
parser.add_argument("-o", "--outdir", dest="outdir", required=True, help="the directory for the output file", metavar = 'OUTDIR')
args = parser.parse_args()

input_file = args.filename
outdir = args.outdir
if not os.path.isdir(outdir):
     os.mkdir(outdir)
        
# Genetic code

code = {'TTT' : 'F', 'TTC' : 'F', 'TTA' : 'L', 'TTG' : 'L', 
        'TCT' : 'S', 'TCC' : 'S', 'TCA' : 'S', 'TCG' : 'S', 
        'TAT' : 'Y', 'TAC' : 'Y', 'TAA' : '*', 'TAG' : '*', 
        'TGT' : 'C', 'TGC' : 'C', 'TGA' : '*', 'TGG' : 'W', 
        'CTT' : 'L', 'CTC' : 'L', 'CTA' : 'L', 'CTG' : 'L', 
        'CCT' : 'P', 'CCC' : 'P', 'CCA' : 'P', 'CCG' : 'P', 
        'CAT' : 'H', 'CAC' : 'H', 'CAA' : 'Q', 'CAG' : 'Q', 
        'CGT' : 'R', 'CGC' : 'R', 'CGA' : 'R', 'CGG' : 'R', 
        'ATT' : 'I', 'ATC' : 'I', 'ATA' : 'I', 'ATG' : 'M', 
        'ACT' : 'T', 'ACC' : 'T', 'ACA' : 'T', 'ACG' : 'T', 
        'AAT' : 'N', 'AAC' : 'N', 'AAA' : 'K', 'AAG' : 'K', 
        'AGT' : 'S', 'AGC' : 'S', 'AGA' : 'R', 'AGG' : 'R', 
        'GTT' : 'V', 'GTC' : 'V', 'GTA' : 'V', 'GTG' : 'V', 
        'GCT' : 'A', 'GCC' : 'A', 'GCA' : 'A', 'GCG' : 'A', 
        'GAT' : 'D', 'GAC' : 'D', 'GAA' : 'E', 'GAG' : 'E', 
        'GGT' : 'G', 'GGC' : 'G', 'GGA' : 'G', 'GGG' : 'G'}

# Reading the input file

def header_remover(file):
    """
    Removes the header from the vep output in tab format, creates a copy of the file without the header
    :param file: path to vep output in tab format
    :returns: the name of the copy of the file without the header
    """
    bashCommand = 'cat ' + file + ' | grep -v ^## > ' + file +'.cl'
    os.system(bashCommand)
    return file + '.cl'

def input_reader(file):
    """
    Reads a file in tab formal without a header
    :param file: path to the input file
    :returns: the data from the file as pandas.DataFrame
    """
    out = pd.read_csv(file, sep='\t', engine='python')
    bashCommand = 'rm ' + file
    os.system(bashCommand) #removing the file without the header
    return out

# Searching for MNVs and annotating them

def close_finder(table):
    """
    Finds variations in the same codon
    :param table: a pandas.DataFrame containing vep annotation data
    :returns: a pandas.DataFrame with variations in the same codon
    """
    out = pd.DataFrame(columns=table.columns)
    table = table.sort_values('Gene')
    genes = list(table['Gene'].unique())
    for i in genes:
        positions = []
        gene_table = table.loc[table['Gene'] == i]
        gene_table = gene_table.sort_values('Protein_position')
        if gene_table.shape[0] >= 2:
            for j in range(gene_table.shape[0]-1):
                if gene_table.iloc[j]['Protein_position'] == gene_table.iloc[j+1]['Protein_position']:
                    if gene_table.iloc[j]['STRAND'] == gene_table.iloc[j+1]['STRAND']:
                        positions.append(gene_table.iloc[j]['Protein_position'])
                        positions.append(gene_table.iloc[j+1]['Protein_position'])
        gene_table = gene_table[gene_table['Protein_position'].isin(positions)]
        out = pd.concat([out, gene_table])
    return out

def MNF_remover(table):
    """
    Removes repeats and misannotated MNVs from the vep annotation data
    :param table: a pandas.DataFrame with vep annotation
    :returns: a pandas.DataFrame with the vep annotation data without repeats and misannotated MNVs
    """
    return pd.concat([table,close_finder(table)]).drop_duplicates(keep=False)

def consequences(table):
    """
    Determines the resulting codon, amino acid, consequences and inpact of the found MNVs
    :param table: a pandas.DataFrame with MNVs
    :returns: annotated pandas.DataFrame with MNVs
    """
    out = pd.DataFrame(columns=table.columns)
    table = table.sort_values('Gene')
    genes = list(table['Gene'].unique())
    
    for i in genes:
        gene_table = table.loc[table['Gene'] == i]
        reference_aa = gene_table.iloc[0]['Amino_acids'][-1]
        reference_codon = list(gene_table.iloc[0]['Codons'].split('/')[1].lower())
        position = gene_table.iloc[0]['Protein_position']
        codon = list(gene_table.iloc[0]['Codons'].split('/')[1].lower())
        for j in range(gene_table.shape[0]):
            mutant_codon = list(gene_table.iloc[j]['Codons'][0:3])
            for k in range(3):
                if mutant_codon[k].isupper():
                    codon[k] = mutant_codon[k]
        score = 0
        for c in range(3):
            if codon[c].isupper() and codon[c].lower() != reference_codon[c].lower():
                score += 1
                
        if score > 1:
            codons = ''.join(codon)+'/'+''.join(reference_codon)
            new_table = gene_table.iloc[:1]
            new_table['Codons'] = codons
            new_aa = code[(''.join(codon)).upper()]
            aa = new_aa + '/' + reference_aa
            new_table['Amino_acid'] = aa
            
            # Calculating consequence
            if new_aa == reference_aa:
                consequence = 'synonymous_variant'
            elif new_aa == '*':
                consequence = 'stop_gained'
            elif reference_aa == '*':
                consequence = 'stop_lost'
            elif position == 1:
                consequence = 'start_lost'
            else:
                consequence = 'missense_variant'
            new_table['Consequence'] = consequence
            
            # Calculating impact
            if consequence == 'synonymous_variant':
                impact = 'LOW'
            elif consequence == 'missense_variant':
                impact = 'MODERATE'
            else:
                impact = 'HIGH'
            new_table['IMPACT'] = impact
            
            out = pd.concat([out, new_table])
        
    return out

# Writing the output file def output tab file

def output_writer(table, path):
    """
    Writes annotation in tab format
    :param table: a pandas.DataFrame with annotation
    :param path: the directory for the output 
    """
    table.to_csv(path, sep ='\t')
    
# Creating a pipeline to combine the functions

def MNV_pipeline(file, outpath):
    """
    Combines reading of the data, searching for polymorphisms in the same codon, annotation of MNVs and saving the results in tab format.
    There is an option to save SNPs annotation without repeats and misannotated MNVs
    :param table: a pandas.DataFrame with annotation
    :param path: the directory for the output
    """
    cleaned_file = header_remover(file)
    print('Reading data...')
    data = input_reader(cleaned_file)
    print('Searching for polymorphisms in the same codon...')
    close_data = close_finder(data)
    snp_out = input('Would like to save file with SNP annotation without misannotated MNVs? y/n ')
    if snp_out == 'y':
        snp_data = MNF_remover(data)
        path = outpath + '/snps.tab'
        output_writer(snp_data, path)
    print('Annotating MNVs...')
    mnv_data = consequences(close_data)
    path = outpath + '/mnvs.tab'
    output_writer(mnv_data, path)
    print('Succesfully finished')
    print('Found', mnv_data.shape[0], 'MNVs:')
    consequence = mnv_data['Consequence'].value_counts().to_dict()
    for key in consequence:
        print(key, consequence[key])
    print('The output file(s) can be found in ', outpath)

# Running the pipeline

warnings.filterwarnings("ignore")
MNV_pipeline(input_file, outdir)