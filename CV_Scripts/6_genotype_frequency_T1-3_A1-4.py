#This program is for calculating frequency of genotypes T, T1-3 and A1-A4
import csv
import sys
import re

FileName = sys.argv[1]
SNPfile=open(FileName,'r')
SNPline=SNPfile.readlines()


newfile="6_"+FileName[0:-6]+"_genotype_frequency.csv"# create a new file for writing data
openfile=open(newfile,"w")

#ref indicates YFP (ancestor) protein sequence
ref="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKLHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*"

#the following codes are used for grouping sequences based on generation of phase II 
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
#the following codes are used for calculating frequency of genotypes T, T1-3 and A1-A4
def replicate_list(numSeq,seqlist):
	freList=[]
	Tall=0
	T=0
	T1=0
	T2=0
	T3=0
	A1=0
	A2=0
	A3=0
	A4=0

	for line in seqlist:
		if line[65]=="S" and line[203]=="C":
			Tall=Tall+1
#calculating frequency of the genotype T2 
			if line[46]=="L"and line[128]!="T" and line[140]!="R":
				T2=T2+1	
#calculating frequency of the genotype T1 			
			if line[64]=="L":
				T1=T1+1
#calculating frequency of the genotype T3 
			if line[42]=="M":
				T3=T3+1
#calculating frequency of the genotype A1 and its derived genotype with one fewer mutation 
			if line[64]=="S"and line[101]=="R" and line[144]=="S"and line[163]=="A":
				A1=A1+1
			if line[64]=="S"and line[144]=="S"and line[163]=="A" and line[101]!="R":
				A1=A1+1
#calculating frequency of the genotype A2 and its derived genotypes with one fewer mutation 
			if line[71]=="I"and line[166]=="E" and line[171]=="V":
				A2=A2+1
			if line[166]=="E" and line[171]=="V" and line[71]!="I":
				A2=A2+1
			if line[71]=="I" and line[166]=="E" and line[171]!="V":
				A2=A2+1
			if line[71]=="I"and line[171]=="V" and line[166]!="E":
				A2=A2+1
#calculating frequency of the genotype A3 
			if line[46]=="L" and line[128]=="T" and line[140]=="R":
				A3=A3+1	
#calculating frequency of the genotype A4 
			if line[71]=="C" and line[167]=="V":
				A4=A4+1
#calculating frequency of the genotype T 
			if line[42]!="M" and line[46]!="L" and line[64]!="S"and line[64]!="L" and line[71]!="I" and line[71]!="C" and line[101]!="R" and line[128]!="T" and line[140]!="R" and line[144]!="S" and line[163]!="A" and line[166]!="E" and line[167]!="V" and line[171]!="V":	
				T=T+1

	if numSeq!=0:
		if T>1:
			F1=T*100.00/numSeq
		else:
			F1=0
		if T1>1:
			F2=T1*100.00/numSeq
		else:
			F2=0
		if T2>1:
			F3=T2*100.00/numSeq
		else:
			F3=0
		if T3>1:
			F4=T3*100.00/numSeq
		else:
			F4=0
		if A1>1:
			F5=A1*100.00/numSeq
		else:
			F5=0
		if A2>1:
			F6=A2*100.00/numSeq
		else:
			F6=0
		if A3>1:
			F7=A3*100.00/numSeq
		else:
			F7=0
		if A4>1:
			F8=A4*100.00/numSeq
		else:
			F8=0
		
		fre1=str("%.3f"%F1)
		fre2=str("%.3f"%F2)
		fre3=str("%.3f"%F3)
		fre4=str("%.3f"%F4)
		fre5=str("%.3f"%F5)
		fre6=str("%.3f"%F6)
		fre7=str("%.3f"%F7)
		fre8=str("%.3f"%F8)

		freList.append(fre1)
		freList.append(fre2)
		freList.append(fre3)
		freList.append(fre4)
		freList.append(fre5)
		freList.append(fre6)
		freList.append(fre7)
		freList.append(fre8)

	return freList

#the following codes are used for writing the result into a csv file
with open(newfile, 'wb') as csvfile:
	spamwriter = csv.writer(csvfile)
	spamwriter.writerow(["T","T1","T2","T3","A1","A2","A3","A4"])

	spamwriter.writerow(replicate_list(n1,seqlist1))
	spamwriter.writerow(replicate_list(n2,seqlist2))
	spamwriter.writerow(replicate_list(n3,seqlist3))
	spamwriter.writerow(replicate_list(n4,seqlist4))
