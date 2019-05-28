#This program is for calculating frequency of the typical genotype T2 and its subsets during evolution
import csv
import sys
import re
import itertools

FileName = sys.argv[1]
SNPfile=open(FileName,'r')
SNPline=SNPfile.readlines()

newfile="13_"+FileName[0:-6]+"_paths_of_genotype_T2.csv"# create a new file for writing data
openfile=open(newfile,"w")

#ref indicates the YFP (ancestor) protein sequence
ref="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKLHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*"

#the following codes are used for grouping sequences based on generation

seqlist1=[]
seqlist2=[]
seqlist3=[]
seqlist4=[]
seqlist5=[]

n1=0
n2=0
n3=0
n4=0
n5=0

for line in SNPline:
	if re.search(">",line):
		seqName=line
	else:
		if re.search("phase-I_4th",seqName):
			seqlist1.append(line.strip())
			n1=n1+1
		elif re.search("phase-II_1st",seqName):
			seqlist2.append(line.strip())
			n2=n2+1
		elif re.search("phase-II_2nd",seqName):
			seqlist3.append(line.strip())
			n3=n3+1
		elif re.search("phase-II_3rd",seqName):
			seqlist4.append(line.strip())
			n4=n4+1
		elif re.search("phase-II_4th",seqName):
			seqlist5.append(line.strip())
			n5=n5+1

#the following codes are used for generating mutation combinations with "-" and "+' respectively indicating an incidence of the corresponding mutation or not
def genotype_generator(muta_num):
	genotype=[]
	snp_combination=list(itertools.product(["-","+"],repeat=muta_num))
	for i in snp_combination:
		genotype.append("".join(i))
	return genotype
#the following codes are used for generating a genotype list
def genotype_generator2(muta_num,pos):
	list_genotype=[]
	for i in itertools.product(["-","+"],repeat=muta_num):
		list_genotype.append(i[pos])
	return list_genotype
#the following codes are used for calculating frequency of each genotype
def replicate_list(numSeq,seqlist):
	freList=[]
	seqref=genotype_generator(3)
	genotype_list=[]
	for line in seqlist:
		seq=""
		if line[65]=="S":
			seq=seq+"+"
		else:
			seq=seq+"-"

		if line[203]=="C":
			seq=seq+"+"
		else:
			seq=seq+"-"

		if line[46]=="L":
			seq=seq+"+"
		else:
			seq=seq+"-"
		genotype_list.append(seq)
	if numSeq!=0:
		for i in range(0,len(seqref)):
			seqnum=genotype_list.count(seqref[i])
			seqfre=seqnum*100.00/numSeq
			freList.append(seqfre)
	return freList

#the following codes are used for writing the result into a csv file
with open(newfile, 'wb') as csvfile:
	spamwriter = csv.writer(csvfile)
	spamwriter.writerow(["G66S"]+genotype_generator2(3,0))
	spamwriter.writerow(["Y204C"]+genotype_generator2(3,1))
	spamwriter.writerow(["F47L"]+genotype_generator2(3,2))
	spamwriter.writerow(["phase-I_4th"]+replicate_list(n1,seqlist1))
	spamwriter.writerow(["phase-II_1st"]+replicate_list(n2,seqlist2))
	spamwriter.writerow(["phase-II_2nd"]+replicate_list(n3,seqlist3))
	spamwriter.writerow(["phase-II_3rd"]+replicate_list(n4,seqlist4))
	spamwriter.writerow(["phase-II_4th"]+replicate_list(n5,seqlist5))
