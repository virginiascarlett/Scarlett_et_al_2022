"""
See what proportion of total counts comes from each subgenome.
Run after running getGeneLengths.py.
Include counts file as a command line arg like so:
python3 analyze_counts_by_subg.py 1_HTScounts.txt
"""
import sys

inF = sys.argv[1]

countsDict = {}  #This will be a dictionary with key=geneID and value=count. 
with open(inF, 'r') as countsF:
    for line in countsF.readlines():
        countsList = line.split('\t')
        if not countsList[0].startswith('__') and not countsList[0].startswith('LT'):   #To accommodate plastid genes and the stats at the end of an HTScounts file, e.g. __no_feature  26675504
            geneID = countsList[0]
            value = int(countsList[1].strip('\n'))
            countsDict[geneID] = value

gene2subg = {} #Dictionary with key=geneID and value=subgenome
#with open('BhybridumBhyb26v2.1.gene_exons.gff3', 'r') as gff:
#For analyzing old ABR113 RNA-seq data:
with open('/global/cscratch1/sd/vstartag/RNA_Seq_2018/refGenomes/Bhybridum_plus_plastid.gff3', 'r') as gff:
	read = gff.readlines()
	#ABR113:
	for n in range(3, len(read)):
	#Bhyb26:
	#for n in range(6, len(read)):
		L = read[n].split('\t')
		if L[2] == 'gene':
			geneID = L[8].split(';')[0].lstrip('ID=')
			subg = L[0][:3]
			gene2subg[geneID] = subg

BhDcounts = 0
BhScounts = 0
for gene, count in countsDict.items():
	if gene2subg[gene] == 'BhD':
		BhDcounts += count
	if gene2subg[gene] == 'BhS':
		BhScounts += count

BhD_total_Tr_lengths = 0
#For analyzing old ABR113 RNA-seq data:
with open('/global/cscratch1/sd/vstartag/RNA_Seq_2018/refGenomes/geneLengths/BhDprimTrlengths.txt', 'r') as BhDlens:
#with open('BhDprimTrlengths.txt', 'r') as BhDlens:
    for line in BhDlens.readlines():
        BhD_total_Tr_lengths += int(line.split('\t')[1].strip('\n'))

BhS_total_Tr_lengths = 0
#For analyzing old ABR113 RNA-seq data:
with open('/global/cscratch1/sd/vstartag/RNA_Seq_2018/refGenomes/geneLengths/BhSprimTrlengths.txt', 'r') as BhSlens:
#with open('BhSprimTrlengths.txt', 'r') as BhSlens:
    for line in BhSlens.readlines():
        BhS_total_Tr_lengths += int(line.split('\t')[1].strip('\n'))

#THIS PART IS MISLEADING. YOU CAN'T COMPARE TRANSCRIPTOME SIZE IN BP VS. COUNTS, BECAUSE COUNTS DON'T ACCOUNT FOR GENE LENGTH
#COMPARE NUMBER OF GENES TO NUMBER OF COUNTS INSTEAD
#print("{}% of the (primary) transcriptome is from BhD transcripts.".format(
#	round(BhD_total_Tr_lengths/(BhS_total_Tr_lengths+BhD_total_Tr_lengths)*100, 3)
#	))
print("{}% of counts are from BhD transcripts.".format(
	round(BhDcounts/(BhDcounts+BhScounts)*100, 3)
	))
print("{} counts from BhD transcripts and {} counts are from BhS transcripts.".format(BhDcounts, BhScounts))
