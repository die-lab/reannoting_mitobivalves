import re
import os
import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import sys

#the script is taken mostrly by the previous version but edited by chat to make it plot different width figures.

# mitochondrial invertebrate start codons from ncbi
alt_codons = ['ATG', 'ATA', 'ATT', 'ATC', 'GTG', 'TTG']
genes_name_list = {
    'ATP6': ['ATP6', 'ATPASE 6', 'ATP SYNTHASE F0 SUBUNIT 6', 'ATP SYNTHASE SUBUNIT 6', 'ATP SYNTHETASE SUBUNIT 6', 'ATPASE SUBUNIT 6'],
    'ATP8': ['ATP SYNTHASE F0 SUBUNIT 8', 'ATP8', 'TRUNCATED ATP SYNTHASE F0 SUBUNIT 8', 'ATP SYNTHASE SUBUNIT 8'],
    'CYTB': ['COB', 'CYTB', 'CYTOCHROME B'],
    'COX1': ['CYTOCHROME C OXIDASE SUBUNIT 1', 'COI', 'COX1', 'COXI', 'CYTOCHROME C SUBUNIT I', 'CYTOCHROME C OXIDASE SUBUNIT I'],
    'COX2': ['CYTOCHROME C OXIDASE SUBUNIT 2', 'COII', 'COX2', 'COXII', 'CYTOCHROME C SUBUNIT II', 'CYTOCHROME OXIDASE C SUBUNIT 2', 'CYTOCHROME C OXIDASE SUBUNIT II'],
    'COX3': ['CYTOCHROME C OXIDASE SUBUNIT III', 'CYTOCHROME C OXIDASE SUBUNIT 3', 'COIII', 'COX3', 'COXIII', 'CYTOCHROME OXIDASE C SUBUNIT 3', 'CYTOCHROME C SUBUNIT III'],
    'FORF': ['F-ORF', 'FORF', 'F-SPECIFIC ORF PROTEIN', 'F ORF'],
    'HORF': ['HORF', 'H OPEN READING FRAME'],
    'MORF': ['M-SPECIFIC MORF PROTEIN', 'M-ORF', 'MORF', 'M-SPECIFIC ORF PROTEIN', 'M-ORF PROTEIN', 'M-ORF1', 'M-ORF2', 'MORF21'],
    'NAD1': ['NAD1', 'ND1', 'NADH DEHYDROGENASE SUBUNIT 1'],
    'NAD2': ['NAD2', 'ND2', 'NADH DEHYDROGENASE SUBUNIT 2'],
    'NAD3': ['NAD3', 'ND3', 'NADH DEHYDROGENASE SUBUNIT 3'],
    'NAD4': ['NAD4', 'ND4', 'NADH DEHYDROGENESE SUBUNIT 4', 'NADH DEHYDROGENASE SUBUNIT 4'],
    'NAD4L': ['NADH DEHYDROGENASE SUBUNIT 4L', 'NADH DEHYDROGENASE SUBUNIT L', 'NAD4L', 'ND4L', 'NDL'],
    'NAD5': ['NAD5', 'ND5', 'NADH DEHYDROGENASE SUBUNIT 5'],
    'NAD6': ['NADH DEHYDROGENASE SUBUNIT 6', 'NAD6', 'ND6']
}
mito_genes = genes_name_list.keys()

# create a dictionary with notes on termination codons and stuff
genbank_file = sys.argv[1]

notes_dict = {}
for gene in genes_name_list.keys():
    notes_dict[gene] = {}
    for record in SeqIO.parse(genbank_file, 'genbank'):
        species_db = dict(zip(list(genes_name_list), [[] for i in list(genes_name_list)]))
        gene_complete_check = []
        for feature in record.features:
            if feature.type == 'CDS':
                gene_name = feature.qualifiers.get('gene', [''])[0].upper()
                gene_name_corr = gene_name
                for i in genes_name_list:
                    if gene_name in genes_name_list[i]:
                        gene_name_corr = i
                    else:
                        gene_name_corr = gene_name
                gene_name = gene_name_corr
                note = feature.qualifiers.get('note', [''])[0].upper()
                if len(note) > 0:
                    notes_dict[gene][record.id] = note

# function to check if length is divisible by 3
def codons_of_three(input):
    if (len(record.seq) / 3 % 1) > 0:
        w = 1
    else:
        w = 0
    return w, len(record.seq) / 3

def split_by_three(string):
    chunks = [string[i:i + 3] for i in range(0, len(string), 3)]
    if len(chunks[-1]) < 3:
        chunks[-2] = chunks[-2] + chunks[-1]
        del chunks[-1]
    return chunks

# function to plot the barplot
def plot_start(alternative_dict, gene, height, weird, notes, no_notes, max_length_parse):
    start_counts = {}
    for codon in alt_codons:
        counts = np.fromiter(alternative_dict[gene][codon].values(), dtype=int)
        positions = np.fromiter(alternative_dict[gene][codon].keys(), dtype=int)
        start_counts[codon] = counts

    fig, ax = plt.subplots(figsize=(max_length_parse * 0.1, 24))  # Adjust width dynamically
    bottom = np.zeros(len(positions))
    for sex, start_count in start_counts.items():
        p = ax.bar(positions, start_count, label=sex, bottom=bottom, width=3)
        bottom += start_count

    ax.set_title(gene)
    ax.axhline(height)
    ax.axhline(height - weird, color='black', ls='--', linewidth=0.5)
    ax.axhline(height - no_notes, color='red', ls='dotted', linewidth=0.5)
    ax.axhline(height - notes, color='green', ls='dotted', linewidth=0.5)
    ax.get_yaxis().set_visible(False)
    plt.ylim((0, numoffiles))
    plt.savefig(f"{gene}_start_codons.svg")

current_directory = os.getcwd()
files_in_directory = os.listdir(current_directory)

# number of species
numoffiles = len(set(''.join(file.split('_')[0:2]) for file in files_in_directory if file.endswith('.fasta')))

# main code
alternative_dict = {}
for gene in mito_genes:
    lookup = f"_{gene}_"
    gene_files = [file for file in files_in_directory if lookup in file and file.endswith('.fasta')]
    if len(gene_files) == 0:
        continue

    weird_length = 0
    no_notes_about = 0
    notes_about = 0
    length_tot = []

    for file in gene_files:
        record = SeqIO.read(file, "fasta")
        cot = codons_of_three(record)
        if cot[0] == 1:
            if record.id in notes_dict[gene]:
                if 'START' in notes_dict[gene][record.id]:
                    notes_about += 1
            else:
                no_notes_about += 1
        length_tot.append(cot[1])
        weird_length += cot[0]

    weird_length_fraction = weird_length / len(gene_files)
    alternative_dict[gene] = {}

    for alternative in alt_codons:
        my_dict = {i: 0 for i in range(0, int(max(length_tot)) * 3)}
        for file in gene_files:
            record = SeqIO.read(file, "fasta")
            positions = [i * 3 for i, chunk in enumerate(split_by_three(str(record.seq))) if chunk == alternative]
            pos_corr = [element for element in positions]
            for x in pos_corr:
                if x < int(max(length_tot)) * 3:
                    my_dict[x] += 1
        alternative_dict[gene][alternative] = my_dict

    plot_start(alternative_dict, gene, len(gene_files), weird_length, notes_about, no_notes_about, max(length_tot))
