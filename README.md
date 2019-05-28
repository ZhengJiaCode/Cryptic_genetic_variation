# Cryptic_genetic_variation

The scripts in the folder "CV_Scripts" were used for SMRT sequencing data analysis in the manuscript titled 'Cryptic genetic variation accelerates evolution by opening access to diverse adaptive peaks'. The SMRT sequencing data are in the folder "SMRT_seq".



The first twenty files are python scripts (Python 2.7.12). Please first install python 2.7.x on your computer to execute these programs. An example for executing programs below (via the terminal of ubuntu 16.04LTS system):

python 1_Mean_numberOfSNPs_cDNA.py v0_replicate1.fasta

Here the input file 'v0_replicate1.fasta' (in fasta format) contains cDNA sequences of YFP individuals sequenced by SMRT sequencing in the population v0 replicate 1 during evolution. Note that you may need specify the paths for '1_Mean_numberOfSNPs_cDNA.py' and 'v0_replicate1.fasta' if they are not in the current directory.



The scripts "21_data2gexf.m" and "22_edge2dom.c" written by Prof. Dr. Joshua Payne were used to produce a .gexf file which we use in Gephi to make an image.
