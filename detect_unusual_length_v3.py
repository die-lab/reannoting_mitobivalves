from Bio import SeqIO
import itertools as it
import os
import numpy as np
from scipy.stats import norm

#instead of printing only genes with length significantly different, print a table with every length, to be analysed in R

#I was trying to select accession numer and genes with significantly different gene length (fist dict with morf, forf and horf for bivalves)
#genes_name_list = { 'ATP6':['ATP6', 'ATPASE 6','ATP SYNTHASE F0 SUBUNIT 6','ATP SYNTHASE SUBUNIT 6','ATP SYNTHETASE SUBUNIT 6','ATPASE SUBUNIT 6'], 'ATP8':['ATP SYNTHASE F0 SUBUNIT 8','ATP8','TRUNCATED ATP SYNTHASE F0 SUBUNIT 8','ATP SYNTHASE SUBUNIT 8'], 'CYTB':['COB','CYTB','CYTOCHROME B'], 'COX1':[',CYTOCHROME C OXIDASE SUBUNIT 1','COI','COX1','COXI','CYTOCHROME C SUBUNIT I','CYTOCHROME C OXIDASE SUBUNIT I'], 'COX2':['CYTOCHROME C OXIDASE SUBUNIT 2','COII','COX2','COXII','CYTOCHROME C SUBUNIT II','CYTOCHROME OXIDASE C SUBUNIT 2','CYTOCHROME C OXIDASE SUBUNIT II'], 'COX3':['CYTOCHROME C OXIDASE SUBUNIT III','CYTOCHROME C OXIDASE SUBUNIT 3','COIII','COX3','COXIII','CYTOCHROME OXIDASE C SUBUNIT 3','CYTOCHROME C SUBUNIT III'], 'FORF':['F-ORF', 'FORF','F-SPECIFIC ORF PROTEIN','F ORF'], 'HORF':['HORF','H OPEN READING FRAME'], 'MORF':['M-SPECIFIC MORF PROTEIN','M-ORF', 'MORF','M-SPECIFIC ORF PROTEIN','M-ORF PROTEIN','M-ORF1','M-ORF2','MORF21'],'NAD1':['NAD1','ND1','NADH DEHYDROGENASE SUBUNIT 1'],'NAD2':['NAD2', 'ND2','NADH DEHYDROGENASE SUBUNIT 2'],'NAD3':['NAD3','ND3','NADH DEHYDROGENASE SUBUNIT 3'],'NAD4':['NAD4', 'ND4','NADH DEHYDROGENESE SUBUNIT 4','NADH DEHYDROGENASE SUBUNIT 4'],'NAD4L':['NADH DEHYDROGENASE SUBUNIT 4L','NADH DEHYDROGENASE SUBUNIT L','NAD4L', 'ND4L','NDL'],'NAD5':['NAD5', 'ND5','NADH DEHYDROGENASE SUBUNIT 5'],'NAD6':['NADH DEHYDROGENASE SUBUNIT 6','NAD6', 'ND6']}
genes_name_list = { 'ATP6':['ATP6', 'ATPASE 6','ATP SYNTHASE F0 SUBUNIT 6','ATP SYNTHASE SUBUNIT 6','ATP SYNTHETASE SUBUNIT 6','ATPASE SUBUNIT 6'], 'ATP8':['ATP SYNTHASE F0 SUBUNIT 8','ATP8','TRUNCATED ATP SYNTHASE F0 SUBUNIT 8','ATP SYNTHASE SUBUNIT 8'], 'CYTB':['COB','CYTB','CYTOCHROME B'], 'COX1':[',CYTOCHROME C OXIDASE SUBUNIT 1','COI','COX1','COXI','CYTOCHROME C SUBUNIT I','CYTOCHROME C OXIDASE SUBUNIT I'], 'COX2':['CYTOCHROME C OXIDASE SUBUNIT 2','COII','COX2','COXII','CYTOCHROME C SUBUNIT II','CYTOCHROME OXIDASE C SUBUNIT 2','CYTOCHROME C OXIDASE SUBUNIT II'], 'COX3':['CYTOCHROME C OXIDASE SUBUNIT III','CYTOCHROME C OXIDASE SUBUNIT 3','COIII','COX3','COXIII','CYTOCHROME OXIDASE C SUBUNIT 3','CYTOCHROME C SUBUNIT III'],'NAD1':['NAD1','ND1','NADH DEHYDROGENASE SUBUNIT 1'],'NAD2':['NAD2', 'ND2','NADH DEHYDROGENASE SUBUNIT 2'],'NAD3':['NAD3','ND3','NADH DEHYDROGENASE SUBUNIT 3'],'NAD4':['NAD4', 'ND4','NADH DEHYDROGENESE SUBUNIT 4','NADH DEHYDROGENASE SUBUNIT 4'],'NAD4L':['NADH DEHYDROGENASE SUBUNIT 4L','NADH DEHYDROGENASE SUBUNIT L','NAD4L', 'ND4L','NDL'],'NAD5':['NAD5', 'ND5','NADH DEHYDROGENASE SUBUNIT 5'],'NAD6':['NADH DEHYDROGENASE SUBUNIT 6','NAD6', 'ND6']}

mito_genes = genes_name_list.keys()

current_directory = os.getcwd()
files_in_directory = os.listdir(current_directory)

p_value = 0.05
def make_distribution(length_list,value):
	cdf_value = norm.cdf(value, loc=np.mean(length_list), scale=np.std(length_list))
	return cdf_value

tab_dict = {}
files = [file.split('_')[0]+"_"+file.split('_')[1] for file in files_in_directory]

for code in set(files):
	tab_dict[code] = {}
	for gene in mito_genes:
		tab_dict[code][gene] = "NA"
		
for gene in mito_genes:
	lookup = str('_'+gene+'_')		
	gene_files = [file for file in files_in_directory if lookup in file and file.endswith('.fasta')]
	#print(gene)
	if len(gene_files) == 0:
		continue
	length_list = []
	for file in gene_files:
		record = SeqIO.read(file,"fasta")
		length_list.append(len(record.seq))
	#distribution_tails = make_distribution(length_list)
	for file in gene_files:
		record = SeqIO.read(file,"fasta")
		count = 0
		tab_dict[record.id][gene] = len(record.seq)

print('mito_code', end='')
for line in mito_genes:
	print(f"\t{line}", end='')
print('/n')

#tab_dict.keys

for line in tab_dict.keys():
	print(line, end='')
	for column in tab_dict[line].values():
		print(f"\t{column}", end = '')
	print()
