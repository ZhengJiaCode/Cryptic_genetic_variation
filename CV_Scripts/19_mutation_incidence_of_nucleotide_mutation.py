#This program is for calculating incidence of each kind of nucleotide mutation
import csv
import sys
import re

SNPfile=open('control_and_ancestor.fasta','r') # open the file which contains the DNA sequences of two 'Control' populations and two 'YFP(ancstor)' populations sequenced by SMRT sequencing
SNPline=SNPfile.readlines()

newfile="19_"+"mutation_incidence_of_nucleotides.csv"# create a new file for writing data
openfile=open(newfile,"w")

#ref indicates the YFP (ancestor) cDNA sequence
ref="ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCTTCGGCTACGGCCTGCAATGCTTCGCCCGCTACCCCGACCACATGAAGCTGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCTACCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAGTGA"

#the following codes are used for grouping sequences from the same replicate population
seqlist1=[]
seqlist2=[]
seqlist3=[]
seqlist4=[]

n1=0
n2=0
n3=0
n4=0


for line in SNPline:
	if re.search(">",line):
		seqName=line
	else:
		if re.search("YFP\(ancestor\)_replicate1",seqName):
			seqlist1.append(line.strip())
			n1=n1+1
		elif re.search("YFP\(ancestor\)_replicate2",seqName):
			seqlist2.append(line.strip())
			n2=n2+1
		elif re.search("Control_replicate1",seqName):
			seqlist3.append(line.strip())
			n3=n3+1
		elif re.search("Control_replicate2",seqName):
			seqlist4.append(line.strip())
			n4=n4+1

#the following codes are used for calculating incidence of each kind of nucleotide mutation
def mutation_list(numSeq,seqlist):
	frequency_list=[]
	A_list=[]
	T_list=[]
	C_list=[]
	G_list=[]
	if numSeq!=0:
		for i in range (0,720):
			for line in seqlist:
				if ref[i]=="A":
					A_list.append(line[i])
				elif ref[i]=="T":
					T_list.append(line[i])
				elif ref[i]=="C":
					C_list.append(line[i])
				elif ref[i]=="G":
					G_list.append(line[i])

		FAA=[A_list.count("A")*1.000000/len(A_list)]
		FAT=[A_list.count("T")*1.000000/len(A_list)]
		FAC=[A_list.count("C")*1.000000/len(A_list)]
		FAG=[A_list.count("G")*1.000000/len(A_list)]

		FTA=[T_list.count("A")*1.000000/len(T_list)]
		FTT=[T_list.count("T")*1.000000/len(T_list)]
		FTC=[T_list.count("C")*1.000000/len(T_list)]
		FTG=[T_list.count("G")*1.000000/len(T_list)]

		FCA=[C_list.count("A")*1.000000/len(C_list)]
		FCT=[C_list.count("T")*1.000000/len(C_list)]
		FCC=[C_list.count("C")*1.000000/len(C_list)]
		FCG=[C_list.count("G")*1.000000/len(C_list)]

		FGA=[G_list.count("A")*1.000000/len(G_list)]
		FGT=[G_list.count("T")*1.000000/len(G_list)]
		FGC=[G_list.count("C")*1.000000/len(G_list)]
		FGG=[G_list.count("G")*1.000000/len(G_list)]
		frequency_list.extend(list(FAA+FAT+FAC+FAG+FTA+FTT+FTC+FTG+FCA+FCT+FCC+FCG+FGA+FGT+FGC+FGG))

	return frequency_list

#the following codes are used for writing the result into a csv file
with open(newfile, 'wb') as csvfile:
	spamwriter = csv.writer(csvfile)
	spamwriter.writerow(["mutation","A->A","A->T","A->C","A->G","T->A","T->T","T->C","T->G","C->A","C->T","C->C","C->G","G->A","G->T","G->C","G->G"])
	spamwriter.writerow(["ancestor_1"]+mutation_list(n1,seqlist1))
	spamwriter.writerow(["ancestor_2"]+mutation_list(n2,seqlist2))
	spamwriter.writerow(["control_1"]+mutation_list(n3,seqlist3))
	spamwriter.writerow(["control_2"]+mutation_list(n4,seqlist4))
