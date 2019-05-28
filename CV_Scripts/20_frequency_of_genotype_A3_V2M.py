#This program is for calculating frequency of genotypes A3, A3+V2M and the ratio between them (r) in the Vc population replicate 3 during the phase II evolution
import csv
import sys
import re

SNPfile=open('2_aa_Vc_replicate3.fasta','r') # open the file '2_aa_Vc_replicate3.fasta' which contains all protein sequences of the Vc population replicate 3 during the phase II evolution
SNPline=SNPfile.readlines()

newfile="20_"+"frequency_of_genotype_A3_V2M.csv"# create a new file for writing data
openfile=open(newfile,"w")

#ref indicates the YFP (ancestor) cDNA sequence
ref="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKLHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*"

#the following codes are used for grouping sequences based on generation
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
		if re.search("phase-II_1st",seqName):
			seqlist1.append(line.strip())
			n1=n1+1
		elif re.search("phase-II_2nd",seqName):
			seqlist2.append(line.strip())
			n2=n2+1
		elif re.search("phase-II_3rd",seqName):
			seqlist3.append(line.strip())
			n3=n3+1
		elif re.search("phase-II_4th",seqName):
			seqlist4.append(line.strip())
			n4=n4+1
#the following codes are used for calculating frequency of genotypes A3, A3+V2M and the ratio between them (r)
def replicate_list(numSeq,seqlist):
	frequency_list=[]
	A3=0
	A3_V2M=0
	for line in seqlist:
		if line[65]=="S" and line[203]=="C":
			if line[46]=="L"and line[128]=="T" and line[140]=="R":
				A3=A3+1
			if line[46]=="L"and line[128]=="T" and line[140]=="R" and line[1]=="M":
				A3_V2M=A3_V2M+1
	if numSeq!=0:
		if A3>0:
			F_A3=A3*100.00/numSeq
		else:
			F_A3=0
		if A3_V2M>0:
			F_A3_V2M=A3_V2M*100.00/numSeq
		else:
			F_A3_V2M=0
		frequency_list.append("%.2f"%F_A3)
		frequency_list.append("%.2f"%F_A3_V2M)
		frequency_list.append("%.2f"%(F_A3_V2M/F_A3*100.0))
	return frequency_list

#the following codes are used for writing the result into a csv file
with open(newfile, 'wb') as csvfile:
	spamwriter = csv.writer(csvfile)
	spamwriter.writerow(["","A3","A3+V2M","r (A3+V2M/A3)"])
	spamwriter.writerow(["phase_II-1"]+replicate_list(n1,seqlist1))
	spamwriter.writerow(["phase_II-2"]+replicate_list(n2,seqlist2))
	spamwriter.writerow(["phase_II-3"]+replicate_list(n3,seqlist3))
	spamwriter.writerow(["phase_II-4"]+replicate_list(n4,seqlist4))

