#This program is for calculating frequency of most populated genotypes during evoloution
import csv
import sys
import re

FileName = sys.argv[1]
SNPfile=open(FileName,'r')
SNPline=SNPfile.readlines()

newfile="11_"+FileName[0:-6]+"_path_of_POM.csv"# create a new file for writing data
openfile=open(newfile,"w")

#ref indicates the YFP (ancestor) protein sequence
ref="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKLHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*"
#the following codes are used for grouping sequences based on generation
seqlist1=[]
seqlist2=[]
seqlist3=[]
seqlist4=[]
seqlist5=[]
seqlist6=[]
seqlist7=[]
seqlist8=[]

n1=0
n2=0
n3=0
n4=0
n5=0
n6=0
n7=0
n8=0

for line in SNPline:
	if re.search(">",line):
		seqName=line
	else:
		if re.search("phase-I_1st",seqName):
			seqlist1.append(line.strip())
			n1=n1+1
		elif re.search("phase-I_2nd",seqName):
			seqlist2.append(line.strip())
			n2=n2+1
		elif re.search("phase-I_3rd",seqName):
			seqlist3.append(line.strip())
			n3=n3+1
		elif re.search("phase-I_4th",seqName):
			seqlist4.append(line.strip())
			n4=n4+1
		elif re.search("phase-II_1st",seqName):
			seqlist5.append(line.strip())
			n5=n5+1
		elif re.search("phase-II_2nd",seqName):
			seqlist6.append(line.strip())
			n6=n6+1
		elif re.search("phase-II_3rd",seqName):
			seqlist7.append(line.strip())
			n7=n7+1
		elif re.search("phase-II_4th",seqName):
			seqlist8.append(line.strip())
			n8=n8+1


#the following codes are used for calculating frequency of  mutants G66S and Y204C, of genotype G66S+Y204C (T), of all high-fitness genotypes that had significantly higher green fluorescence than genotype T, as well as of all intermediate genotypes (averaged) leading to these high fitness genotypes that are inaccessible through selection for green fluorescence alone.
def replicate_list(numSeq,seqlist):
	freList1=[]
	freList2=[]
	freList3=[]
	freList4=[]

	num_intermediate=0 
	num_high_fitness=0
	num_g66s_y204c=0
	num_g66s_or_y204c=0
	for line in seqlist:
		if line[65]=="S" and line[203]=="C":
			if line[46]=="L":
				num_high_fitness=num_high_fitness+1	
			elif line[42]=="M":
				num_high_fitness=num_high_fitness+1				
			elif line[64]=="L":
				num_high_fitness=num_high_fitness+1
			elif line[71]=="C":
				num_high_fitness=num_high_fitness+1
			elif line[71]=="I":
				num_high_fitness=num_high_fitness+1				
			elif line[171]=="V":
				num_high_fitness=num_high_fitness+1
			elif line[64]=="S" and line[163]=="A" and line[144]=="S":
				num_high_fitness=num_high_fitness+1
			elif line[46]!="L" and line[64]!="S" and line[71]!="C" and line[71]!="I" and line[101]!="R" and line[128]!="T" and line[140]!="R" and line[163]!="A" and line[171]!="V" and line[64]!="L" and line[144]!="S" and line[42]!="M" and line[166]!="E" and line[167]!="V":
				num_g66s_y204c = num_g66s_y204c+1
				
		else:
			if line[46]=="L" or line[64]=="S" or line[71]=="C" or line[71]=="I" or line[101]=="R" or line[128]=="T" or line[140]=="R" or line[163]=="A" or line[171]=="V" or line[42]=="M"or line[166]=="E" or line[144]=="S" or line[167]=="V":
				num_intermediate=num_intermediate+1
			else:
				if line[65]=="S" or line[203]=="C":
					num_g66s_or_y204c=num_g66s_or_y204c+1

	if numSeq!=0:
		F1=num_high_fitness*100.00/numSeq
		F2=num_intermediate*100.00/numSeq
		F3=num_g66s_y204c*100.00/numSeq
		F4=num_g66s_or_y204c*100.00/numSeq

		fre1=str("%.3f"%F1)
		fre2=str("%.3f"%F2)
		fre3=str("%.3f"%F3)
		fre4=str("%.3f"%F4)

		freList1.append(fre1)
		freList1.append(fre2)
		freList1.append(fre3)
		freList1.append(fre4)

	return freList1
	return freList2
	return freList3
	return freList4

#the following codes are used for writing the results into a csv file
with open(newfile, 'wb') as csvfile:
	spamwriter = csv.writer(csvfile)
	spamwriter.writerow(["High-fitness genotypes","Intermediate inaccessible genotypes","T","G66S or Y204C"])
	spamwriter.writerow(replicate_list(n1,seqlist1))
	spamwriter.writerow(replicate_list(n2,seqlist2))
	spamwriter.writerow(replicate_list(n3,seqlist3))
	spamwriter.writerow(replicate_list(n4,seqlist4))
	spamwriter.writerow(replicate_list(n5,seqlist5))
	spamwriter.writerow(replicate_list(n6,seqlist6))
	spamwriter.writerow(replicate_list(n7,seqlist7))
	spamwriter.writerow(replicate_list(n8,seqlist8))

