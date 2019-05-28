#This program is for determining distribution of the number of amino acid changing mutations per individual molecule in V0 and VC populations during the phase II evolution

import csv
import sys
import re
FileName = sys.argv[1]
SNPfile=open(FileName,'r')
SNPline=SNPfile.readlines()

newfile="3_"+FileName[0:-4]+"_aa_snp_distribution.csv" # create a new file for reading data
openfile=open(newfile,"w")

#ref indicates YFP (ancestor) protein sequence
ref="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKLHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*"

seqlist1=[]
seqlist2=[]
seqlist3=[]
seqlist4=[]

n1=0
n2=0
n3=0
n4=0

#the following codes are used for grouping sequences based on generation 
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

#the following codes are used for calculating mean number of amino-acid changing mutations per YFP molecule in each replicate population and in each generation 
def replicate_list(numSeq,seqlist):
	freList=[]	
	diffList=[]
	for line in seqlist:
		diff=0
		for i in range (0,240):
			if line[i]!=ref[i]:
				diff=diff+1
		diffList.append(diff)
	for j in range (0, 50):#here we count the number of variants which have amino-acid changing mutations varying from 0 to 50 (note that no variant has more than 50 amino-acid changing mutations).
		Num=diffList.count(j)
		if numSeq!=0:
			F=Num*100.00/numSeq
			freList.append(F)
	return freList

#the following codes are used for writing the result into a csv file
with open(newfile, 'wb') as csvfile:
	spamwriter = csv.writer(csvfile)
	spamwriter.writerow(list(range(0,50)))
	spamwriter.writerow(replicate_list(n1,seqlist1))
	spamwriter.writerow(replicate_list(n2,seqlist2))
	spamwriter.writerow(replicate_list(n3,seqlist3))
	spamwriter.writerow(replicate_list(n4,seqlist4))
