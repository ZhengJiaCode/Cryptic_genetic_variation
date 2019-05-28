#This program is used for calculating frequency of each amino-acid changing mutation that reached more than 10% in at least one replicate population during the phase II evolution
import csv
import sys
import re

FileName = sys.argv[1]
SNPfile=open(FileName,'r')
SNPline=SNPfile.readlines()

newfile="8_"+FileName[0:-6]+"_dynamic_of_aaSNP_frequency>10%.csv"# create a new file for writing data
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

#the following codes are used for calculating frequency of the amino-acid changing mutations that reached >10% in at least one replicte population
replicate_voc={n1:"phase-II_1st",n2:"phase-II_2nd",n3:"phase-II_3rd",n4:"phase-II_4th"}

mutaConvert={1:"A",2:"N",4:"S",6:"K",10:"R",19:"N",42:"M",46:"L",59:"S",64:"L",65:"S",71:"C",73:"H",90:"A",99:"S",101:"R",128:"T",140:"R",
144:"S",152:"V",158:"M",163:"A",164:"S",166:"E",167:"V",171:"T",187:"S",190:"G",192:"S",198:"D",203:"C",223:"L",237:"C"}
mutaConvert2={1:"M",64:"S",71:"I",171:"V",198:"S"}

def replicate_list(numSeq,seqlist):
	if numSeq!=0:
		freList=[]
		freList.append(replicate_voc[numSeq])
		for i in [1,2,4,6,10,19,42,46,59,64,65,71,73,90,99,101,128,140,144,152,158,163,164,166,167,171,187,190,192,198,203,223,237]:
			snpMatrix=[]
			for line in seqlist:
				snpMatrix.append(line[i])
			x=mutaConvert[i]
			Num=snpMatrix.count(x)
			F=Num*100.00/numSeq
			aaSNP=str("%.3f"%F)
			freList.append(aaSNP)

			if i in [1,64,71,171,198]:
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
	for i in [1,2,4,6,10,19,42,46,59,64,65,71,73,90,99,101,128,140,144,152,158,163,164,166,167,171,187,190,192,198,203,223,237]:
		x=mutaConvert[i]
		aaSNP=ref[i]+str(i+1)+x
		snpList.append(aaSNP)
		if i in [1,64,71,171,198]:
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
