# reannoting_mitobivalves
Methods and plots of the project!

This repository wants to be a step-by-step guide of the analyses, while giving you the chance to run it again by your own, custom data.
I suggest to `git clone` the current repository, 'cause you will need some of its scritps.

## download and build datasets
I have searched for mitochondrial genomes on ncbi by command line, using the below command.
You can use different classes of organisms, just editing the *metazoa_list.txt*.
It would be a good idea to see how many mitochondrial we are going to downlaod checking that type of query on the browser [browser](https://www.ncbi.nlm.nih.gov/).
```
for class in $(cat metazoa_list.txt); do esearch -db nuccore -query "("$class"[Organism] AND mitochondrion[All Fields] AND complete[All Fields]" | efetch -format gb >> $class'.gb'; done 
```
For the pourpose of our project, I chose to shrink the number of genes down to only RefSeq sequence. In this way, I wanted to limit stats only on mitochondrial genomes teorically well annotated.

## extract single genes
I have extracted every gene's sequence (nucleotide and protein type) from gb files using python script *extracting_feature_v6.py*. I recommend to keep the next steps into separate directories.
```
for gb_file in $(ls *.gb); do python extracting_feature_v6.py $gb_file; mkdir ${gb_file%.gb}_extraction; mv *.fasta ${gb_file%.gb}_extraction; done
```
Annotations do not follow a golden rule, so I have used a dictionary to assign every name to the 13 mitochondrial genes. You can see it a the begininnig *extracting_feature_v6.py*. Your data can have some genes named differently. I suggest to check manually for bad annotation.

For each group of organisms, move the two type of sequences to different places.
```
for i in *_extraction; do cd $i; mkdir $i'_nucl'; mkdir $i'_prot'; mv *_nucl.fasta $i'_nucl'; mv *_prot.fasta $i'_prot'; cd ..; done
```









