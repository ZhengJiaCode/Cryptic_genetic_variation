#This program is used for getting all genotypes consisted of residues at the below fifteen positions and their frequency, and the generated data are used for creating genotype networks using Gephi
import csv
import sys
import re

FileName = sys.argv[1]
SNPfile=open(FileName,'r')
SNPline=SNPfile.readlines()

newfile_txt="10_"+FileName[0:-4]+"_genotype_network.txt"# create a new file for writing data
openfile_txt=open(newfile_txt,"w")

#ref indicates the YFP (ancestor) protein sequence, and snp10ref indicates the ancestral genotype consisted of the amino acid residues at the following fifteen positions in the ancestor protein sequence
ref="MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFGYGLQCFARYPDHMKLHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSYQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYK*"
snp10ref="VLFFGFKIKNVKIIY"

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
		seq=""
		for i in["2","43","47","65","66","72","102","129","141","145","164","167","168","172","204"]:
			seq=seq+line[int(i)-1]
		if not "*" in seq:
			if re.search("phase-I_1st",seqName):
				seqlist1.append(seq)
				n1=n1+1
			elif re.search("phase-I_2nd",seqName):
				seqlist2.append(seq)
				n2=n2+1
			elif re.search("phase-I_3rd",seqName):
				seqlist3.append(seq)
				n3=n3+1
			elif re.search("phase-I_4th",seqName):
				seqlist4.append(seq)
				n4=n4+1
			elif re.search("phase-II_1st",seqName):
				seqlist5.append(seq)
				n5=n5+1
			elif re.search("phase-II_2nd",seqName):
				seqlist6.append(seq)
				n6=n6+1
			elif re.search("phase-II_3rd",seqName):
				seqlist7.append(seq)
				n7=n7+1
			elif re.search("phase-II_4th",seqName):
				seqlist8.append(seq)
				n8=n8+1


position_index={"0":"2","1":"43","2":"47","3":"65","4":"66","5":"72","6":"102","7":"129","8":"141","9":"145","10":"164","11":"167","12":"168","13":"172","14":"204"}


#the following codes are used for getting the list of all genotypes consisted of residues at the above fifteen positions and their frequency
def genotype_freq(number_Of_Seq,seqlist,group):
	if number_Of_Seq!=0:
		sequence_list=[]
		for line in seqlist:
			sequence_list.append(line)
		unique_sequence=set(sequence_list)
		for x in unique_sequence:
			if not "*" in x:
				Num=sequence_list.count(x)
				Frequency_of_genotype=Num*100.00/number_Of_Seq
				genotypelist=group+"	"+x+"	"+str(Frequency_of_genotype)+"	"+"0"
				openfile_txt.write("%s\n"%genotypelist)


#the following codes are used for writing the results into a txt file which can be used for creating genotype networks using Gephi
openfile_txt.write("Genotypeset	Genotype	Score	Delta\n")
openfile_txt.write("\n")
openfile_txt.write("0	VLFFGFKIKNVKIIY	100.0	0\n")
genotype_freq(n1,seqlist1,"phase-I_1st")
genotype_freq(n2,seqlist2,"phase-I_2nd")
genotype_freq(n3,seqlist3,"phase-I_3rd")
genotype_freq(n4,seqlist4,"phase-I_4th")

genotype_freq(n5,seqlist5,"phase-II_1st")
genotype_freq(n6,seqlist6,"phase-II_2nd")
genotype_freq(n7,seqlist7,"phase-II_3rd")
genotype_freq(n8,seqlist8,"phase-II_4th")
