# Cryptic_genetic_variation

The first twenty python scripts (Python 2.7.12) are used for SMRT sequencing data analysis in the manuscript titled 'Cryptic genetic variation accelerates evolution by opening access to diverse adaptive peaks'.

Please first install python 2.7.x on your computer.
An example for executing programs below (via the terminal of ubuntu 16.04LTS system.Note that you may need specify the paths for '1_Mean_numberOfSNPs_cDNA.py' and 'v0_replicate1.txt' if they are not in the current directory.):

python 1_Mean_numberOfSNPs_cDNA.py v0_replicate1.txt

Here the input file 'v0_replicate1.txt' (in fasta format) contains cDNA sequences of YFP individuals sequenced by SMRT sequencing in the population v0 replicate 1 during evolution. All SMRT sequencing data for this project has been deposited at DDBJ/EMBL/GenBank under the accession KCZY00000000. 
