#This program is used for calculating frequency of 17 adaptive amino-acid changing mutations during the phase I evolution
import csv
import sys
import re

FileName = sys.argv[1]
SNPfile=open(FileName,'r')
SNPline=SNPfile.readlines()


newfile="9_"+FileName[0:-4]+"_dynamic_of_17aaSNP_frequency_phaseI.csv"# create a new file for writing data
openfile=open(newfile,"w")

#ref indicates YFP (ancestor) protein sequence
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

#the following codes are used for calculating frequency of 17 adaptive amino-acid changing mutations in a population in each generation of phase I
replicate_voc={n1:"phase-I_1st",n2:"phase-I_2nd",n3:"phase-I_3rd",n4:"phase-I_4th"}

mutaConvert={1:"A",46:"L",64:"L",65:"S",71:"C",101:"R",128:"T",140:"R",
144:"S",163:"A",166:"E",167:"V",171:"V",203:"C"}
mutaConvert2={1:"M",64:"S",71:"I"}

def replicate_list(numSeq,seqlist):
	if numSeq!=0:
		freList=[]
		freList.append(replicate_voc[numSeq])
		for i in [1,46,64,65,71,101,128,140,144,163,166,167,171,203]:
			snpMatrix=[]
			for line in seqlist:
				snpMatrix.append(line[i])
			x=mutaConvert[i]
			Num=snpMatrix.count(x)
			F=Num*100.00/numSeq
			aaSNP=str("%.3f"%F)
			freList.append(aaSNP)

			if i in [1,64,71]:
				y=mutaConvert2[i]
				Num=snpMatrix.count(y)
				F=Num*100.00/numSeq
				aaSNP=str("%.3f"%F)
				freList.append(aaSNP)

	else:
		freList="0"

	return freList

# the following codes are used for getting the list of amino-acid mutations
def snp_list():
	snpList=["SNP"]
	for i in [1,46,64,65,71,101,128,140,144,163,166,167,171,203]:
		x=mutaConvert[i]
		aaSNP=ref[i]+str(i+1)+x
		snpList.append(aaSNP)
		if i in [1,64,71]:
			y=mutaConvert2[i]
			aaSNP=ref[i]+str(i+1)+y
			snpList.append(aaSNP)
	return snpList
#the following codes are used for writing the results into a csv file
with open(newfile, 'wb') as csvfile:
	spamwriter = csv.writer(csvfile)

	spamwriter.writerow(snp_list())
	spamwriter.writerow(replicate_list(n1,seqlist1))
	spamwriter.writerow(replicate_list(n2,seqlist2))
	spamwriter.writerow(replicate_list(n3,seqlist3))
	spamwriter.writerow(replicate_list(n4,seqlist4))
