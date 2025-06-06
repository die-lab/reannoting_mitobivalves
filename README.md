![Image](https://github.com/user-attachments/assets/ab5ffcfa-1827-4d46-bcfc-19a1ed488c24)

# reannoting_mitobivalves
Methods of the project!

This repository wants to be a step-by-step guide of the analyses, while giving you the chance to run it again by your own, custom data.
I suggest to `git clone` the current repository, 'cause you will need some of its scritps.

## Download and build datasets
I have searched for mitochondrial genomes on ncbi by command line, using the below command.
You can use different classes of organisms, just editing the *metazoa_list.txt*.
It would be a good idea to see how many mitochondrial we are going to downlaod checking that type of query on the [browser](https://www.ncbi.nlm.nih.gov/).
```
for class in $(cat metazoa_list.txt); do esearch -db nuccore -query "("$class"[Organism] AND mitochondrion[All Fields] AND complete[All Fields]" | efetch -format gb >> $class'.gb'; done 
```
For the pourpose of our project, I chose to shrink the number of genes down to only RefSeq sequence. In this way, I wanted to limit stats only on mitochondrial genomes teorically well annotated.

## Extract single genes
I have extracted every gene's sequence (nucleotide and protein type) from gb files using python script *extracting_feature_v6.py*. I recommend to keep the next steps into separate directories.
```
for gb_file in $(ls *.gb); do python extracting_feature_v6.py $gb_file; mkdir ${gb_file%.gb}_extraction; mv *.fasta ${gb_file%.gb}_extraction; done
```
Annotations do not follow a golden rule, so I have used a dictionary to assign every name to the 13 mitochondrial genes. You can see it a the begininnig of *extracting_feature_v6.py*. Your data can have some genes named differently. I also suggest to check manually for bad annotation. Although being uncommon, some mitochondrial genomes can have multiple copy of a genes too. It happend for a couple of genes in bivalves and for complex III,IV and IV of Oegopsida, in cephalopods. While It has been easy to identify the "real" copy in bivalves, in cephalopods they are so conserved that it would have been impossible to find the working one (even because they are both expressed), so I have considered only the first met on annotation. 

For each group of organisms, move the two type of sequences to different places.
```
for i in *_extraction; do cd $i; mkdir $i'_nucl'; mkdir $i'_prot'; mv *_nucl.fasta $i'_nucl'; mv *_prot.fasta $i'_prot'; cd ..; done
```

## Start and stop codons plotting
I wanted to see what was the distribution of start codons along genes, taking the in frame reading. I have used *start_codons_check_v5* and *stop_codons_check_v2*. You can run everything together using the below command. You have to run it in the home directory of the project, where you have *[class].gb* files, *[class]_extraction directories* and scripts. 
```
for i in *extraction; do cd $i; mkdir $i'_start_codons_distribution'; mkdir $i'_stop_codons_distribution'; cd $i'_nucl'; python ../../start_codons_check_v5.py ../../${i%_extraction}.gb > ../$i'_start_codons_distribution'/$i'_start.log'; mv *.svg ../$i'_start_codons_distribution'/.;python ../../stop_codons_check_v2.py ../../${i%_extraction}.gb > ../$i'_stop_codons_distribution'/$i'_stop.log'; mv *.svg ../$i'_stop_codons_distribution'/.; cd ../../.; done
```

## Length and start codons datasets
Get every genes length and start codon usage of your mitochondrial datasets using *detect_unusual_length_v3.py* and *detect_start_codon_usage.py*.
```
for i in *extraction; do cd $i; mkdir $i'_start_codons_distribution'; mkdir $i'_stop_codons_distribution'; cd $i'_nucl'; python ../../detect_unusual_length_v3.py > ../${i%_extraction}_length_dataset.tsv; python ../../detect_start_codon_usage.py > ../${i%_extraction}_start_codon_dataset.tsv; cd ../../.; done
```

## Tree topology and reordering
I have run ML tree inference with modelfinder and 100 bootstraps on align, trim and concatenated genes for each group. I need it mostly because I wanted to see if genes' length and start codons were phylogenetically conserved, at least at the order level. Once You have the tree you have to retrieve the list of species as displayed plotting the tree. You need to reorder each *[class]_length_datasets.tsv* and *[class]_start_codon_dataset.tsv* accordingly to its tree. Now plot everything together!













