#This program is for getting mean number of amino-acid changing mutations per YFP molecule in each replicate population and in each generation
import csv
import sys
import re

FileName = sys.argv[1]
SNPfile=open(FileName,'r')
SNPline=SNPfile.readlines()

newfile="4_"+FileName[0:-6]+"_mean_number_SNPs_aa2.csv" # create a new file for reading data
openfile=open(newfile,"w")
#ref indicates YFP (ancestor) protein sequence
ref="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKLHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*"

seqlist1=[]
seqlist2=[]
seqlist3=[]
seqlist4=[]
seqlist5=[]
seqlist6=[]
seqlist7=[]
seqlist8=[]
seqlist9=[]
seqlist10=[]
seqlist11=[]
seqlist12=[]
n1=0
n2=0
n3=0
n4=0
n5=0
n6=0
n7=0
n8=0
n9=0
n10=0
n11=0
n12=0

#the following codes are used for grouping sequences based on generation 
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
		elif re.search("YFP\(ancestor\)_replicate1",seqName):
			seqlist9.append(line.strip())
			n9=n9+1
		elif re.search("YFP\(ancestor\)_replicate2",seqName):
			seqlist10.append(line.strip())
			n10=n10+1
		elif re.search("Control_replicate1",seqName):
			seqlist11.append(line.strip())
			n11=n11+1
		elif re.search("Control_replicate2",seqName):
			seqlist12.append(line.strip())
			n12=n12+1

#the following codes are used for getting mean number of amino-acid changing mutations per YFP molecule in each replicate population and in each generation 
def replicate_list(numSeq,seqlist):
	freList=[]
	snpNum=0.00
	if numSeq!=0:
		for i in range (0,240):
			diff=0
			for line in seqlist:
				if line[i]!=ref[i]:
					diff=diff+1
			snpNum=snpNum+diff*1.00/numSeq
	freList.append(snpNum)
	return freList

#the following codes are used for writing the result into a csv file
with open(newfile, 'wb') as csvfile:
	spamwriter = csv.writer(csvfile)
	spamwriter.writerow(replicate_list(n1,seqlist1))
	spamwriter.writerow(replicate_list(n2,seqlist2))
	spamwriter.writerow(replicate_list(n3,seqlist3))
	spamwriter.writerow(replicate_list(n4,seqlist4))
	spamwriter.writerow(replicate_list(n5,seqlist5))
	spamwriter.writerow(replicate_list(n6,seqlist6))
	spamwriter.writerow(replicate_list(n7,seqlist7))
	spamwriter.writerow(replicate_list(n8,seqlist8))
	spamwriter.writerow(replicate_list(n9,seqlist9))
	spamwriter.writerow(replicate_list(n10,seqlist10))
	spamwriter.writerow(replicate_list(n11,seqlist11))
	spamwriter.writerow(replicate_list(n12,seqlist12))
