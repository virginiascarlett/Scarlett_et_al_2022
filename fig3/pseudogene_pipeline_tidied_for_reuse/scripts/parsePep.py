"""
Parses peptide sequences from Phytozome from one big file into many little files.
Only produces little files for genes in missing_gene_intervals.sorted.bed.

Include the progenitor(s) on the command line. E.g.:
python3 parsePep.py Bdistachyon_556_v3.2.protein_primaryTranscriptOnly.fa Bstacei_316_v1.1.protein_primaryTranscriptOnly.fa

NEW: One progenitor is also acceptable, e.g.:
python3 parsePep.py Bstacei_316_v1.1.protein_primaryTranscriptOnly.fa

"""

import sys
import glob

def read_fasta(filename): 
	lines = {}
	with open(filename, 'r') as fastafile:
		currentKey = ''
		for line in fastafile.readlines():
			if line.startswith('>'):
				currentKey = line.split()[3][6:] #get locus ID
			else:
				if currentKey in lines:
					lines[currentKey] += line
				else:
					lines[currentKey] = line
	return(lines)

def parse_conserved_gene_info():
	mydata = {}
	with open('conservedgenes.txt', 'r') as myfile:
		read = myfile.readlines()
		for n in range(1, len(read)):
			L = read[n].split('\t')
			if len(L) == 4: #two progenitors
				PP1, PP2, dipgene1, dipgene2 = L[0], L[1], L[2], L[3].strip('\n')
				mydata[PP1] = dipgene1
				mydata[PP2] = dipgene2
			elif len(L) == 3: #one progenitor
				PP1, PP2, dipgene = L[0], L[1], L[2].strip('\n')
				mydata[PP1] = dipgene
				mydata[PP2] = dipgene
	return(mydata)


pepfiles = sys.argv[1:]

pepdata = [] #a list with one dictionary per progenitor
for pepfile in pepfiles:
	pepdata.append(read_fasta(pepfile))

allpeps = pepdata[0] # a single dictionary with all loci from all progenitors
#keys are progenitor gene IDs and values are peptide sequences
# as usual, we are assuming two progenitors maximum
if len(pepdata) > 1:
	allpeps.update(pepdata[1])

bedfilelist = glob.glob('*_gene_intervals.sorted.bed') #can be missing_gene_intervals.sorted.bed or conserved_gene_intervals.sorted.bed
if len(bedfilelist) > 1:
    raise Exception("You must have only one file ending with _gene_intervals.sorted.bed in your working directory.")

bedfilename = bedfilelist[0]

GOIs = [] # (diploid) genes of interest
if bedfilename == 'missing_gene_intervals.sorted.bed':
	with open(bedfilelist[0], 'r') as bedfile:
		for line in bedfile.readlines():
			GOIs.append(line.split("\t")[3].strip("\n"))
elif bedfilename == 'conserved_gene_intervals.sorted.bed':
	PPtodiploid = parse_conserved_gene_info()
	for dipgene in set(PPtodiploid.values()):
		GOIs.append(dipgene)

for GOI in GOIs:
	outF = open('candidatePseudogeneRegions/{}.pep.fa'.format(GOI), 'w')
	outF.write('>{}\n'.format(GOI))
	outF.write(allpeps[GOI])
