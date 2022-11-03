"""
Obtains CDS sequences for diploid genes of interest from phytozome CDS file(s).

Run like so:
python3 parseCDS.py Bdistachyon_556_v3.2.cds_primaryTranscriptOnly.fa Bstacei_316_v1.1.cds_primaryTranscriptOnly.fa

One progenitor is also acceptable, e.g.:
python3 parseCDS.py Bstacei_316_v1.1.cds_primaryTranscriptOnly.fa

"""

import sys
import os
from collections import defaultdict

def read_fasta(filename): 
	lines = {}
	with open(filename, 'r') as fastafile:
		currentKey = ''
		for line in fastafile.readlines():
			if line.startswith('>'):
				currentKey = line.split()[0].strip('>').split(".")[0] #scaffolds won't have unique keys but whatever dgaf
			else:
				if currentKey in lines:
					lines[currentKey] += line.strip('\n')
				else:
					lines[currentKey] = line.strip('\n')
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


cdsfilenames = sys.argv[1:]

superCDSdict = {} #one dict with all the CDS gene IDs as keys and sequences as values
CDSdicts = []
for f in cdsfilenames:
	CDSdicts.append(read_fasta(f))
for d in CDSdicts:
	for k, v in d.items():
		superCDSdict[k] = v

GOIs = [] # (diploid) genes of interest
with open('length_differences.txt', 'r') as myfile:
	for line in myfile.readlines():
		GOIs.append(line.split('\t')[0])

#We need to know which diploid gene IDs correspond to which aliases
#so that each CDS file name will have the appropriate alias.

aliases = {} # If we are working with missing genes, the keys are diploid genes.
#If we are working with conserved genes, the keys are polyploid genes. Values are their aliases
with open('aliases.txt', 'r') as aliasfile:
	for line in aliasfile.readlines():
		alias, geneID = line.split('\t')[0], line.split('\t')[1].strip('\n')
		aliases[geneID] = alias

PPtodiploid = {}
if os.path.exists('conservedgenes.txt'):
	PPtodiploid = parse_conserved_gene_info()

#conserved genes only
if PPtodiploid: #if we are working with conserved genes:
	dipgenes2aliases = defaultdict(list) #keys are diploid gene IDs, value is a list with 
#one alias (if two progenitors) or two aliases (if one progenitor)
	for PPgene, alias in aliases.items():
		mydipgene = PPtodiploid[PPgene]
		dipgenes2aliases[mydipgene].append(alias)
	for dipgene, aliases in dipgenes2aliases.items():
		for a in aliases:
			with open('candidatePseudogeneRegions/{}.cds.fa'.format(a), 'w') as outF:
				outF.write('>{}\n'.format(dipgene))
				outF.write('{}\n'.format(superCDSdict[dipgene]))
	quit() # Stop running the script. We're done.

# missing genes only
for geneID in GOIs:
	with open('candidatePseudogeneRegions/{}.cds.fa'.format(aliases[geneID]), 'w') as outF:
		outF.write('>{}\n'.format(geneID))
		outF.write('{}\n'.format(superCDSdict[geneID]))