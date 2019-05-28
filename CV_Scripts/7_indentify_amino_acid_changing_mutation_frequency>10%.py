#This program is for identifying those amino-acid changing mutations with frequency more than 10% during the phase II evolution
import csv
import sys
import re

FileName = sys.argv[1]
SNPfile=open(FileName,'r')
SNPline=SNPfile.readlines()

newfile="7_"+FileName[0:-4]+"_aaSNP_frequency>10%.csv"# create a new file for writing data
openfile=open(newfile,"w")

#ref indicates YFP (ancestor) protein sequence
ref="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKLHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*"

#the following codes are used for grouping sequences from each generation of phase II evolution 
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

#calculating the frequency of the amino-acid changing mutations with frequency more than 10%
def replicate_list(numSeq,seqlist):
	freList=[]
	for i in range (0,240):
		snpMatrix=[]
		for line in seqlist:
			snpMatrix.append(line[i])
		SetSNP=set(snpMatrix)
		snpList=[]
		for x in SetSNP:
			snpList.append(x)
		for j in range (0, len(snpList)):
			y=snpList[j]
			if y!=ref[i]:
				Num=snpMatrix.count(y)
				F=Num*100.00/numSeq
				if F>=10:
					aaSNP=str(ref[i])+str(i+1)+str(y)+str("%.3f"%F)
					freList.append(aaSNP)
	return freList

#the following codes are used for writing the results into a csv file
with open(newfile, 'wb') as csvfile:
	spamwriter = csv.writer(csvfile)
	spamwriter.writerow(replicate_list(n1,seqlist1))
	spamwriter.writerow(replicate_list(n2,seqlist2))
	spamwriter.writerow(replicate_list(n3,seqlist3))
	spamwriter.writerow(replicate_list(n4,seqlist4))
