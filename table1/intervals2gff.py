"""
Create a new gff file that includes the pseudogenes, as if they were real annotated genes.
We will use this to get counts and TPMs from the RNA-seq libraries.
You can quickly see the pseudogenes in the gff by looking at column 2,
which says "pseudogene" instead of "JGI".
The pseudogene ID is just the geneID for the orthologous diploid gene.
"""

import re

def natural_sort(l): 
	convert = lambda text: int(text) if text.isdigit() else text.lower()
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	return sorted(l, key=alphanum_key)

geneInfo = {} #key is geneID of diploid ortholog, value is a list of lists. 
#each sublist contains coordinates for one exon of the Bhyb26 pseudogene
with open('/global/cscratch1/sd/vstartag/Bhyb26_lost_genesNEWNEW/lost_genes/candidatePseudogeneRegions/pseudogeneIntervals.txt', 'r') as inF:
	for line in inF.readlines():
		L = line.strip('\n').split('\t')
		dipgene, chrom, start, end, strand = L[0], L[5], int(L[6]), int(L[7]), L[8]
		if dipgene in geneInfo:
			geneInfo[dipgene].append([chrom, start, end, strand])
		else:
			geneInfo[dipgene] = [[chrom, start, end, strand]]

results = {} #lines to write to gff file. key looks like e.g. 'BhS1_6311942' where that number is the gene start
for dipgene, masterlist in geneInfo.items():
	newmaster = sorted(masterlist, key=lambda sublist: sublist[1]) #I think the intervals should be sorted by start coord already, but just in case they aren't...
	#it looks like CDSs on the minus strand are written in descending order in the original gff... whatever I'm
	#just going to write all of mine in ascending order and I think it will probably be fine
	chrom = newmaster[0][0]
	genestart = newmaster[0][1]
	geneend = newmaster[-1][2]
	keyID = '{}_{}'.format(chrom, str(genestart))
	results[keyID] = []
	results[keyID].append('{}\tpseudogene\tgene\t{}\t{}\t.\t{}\t.\tID={};Name={}\n'.format(newmaster[0][0], str(genestart), str(geneend), newmaster[0][3], dipgene, dipgene))
	counter = 1
	for sublist in newmaster:
		results[keyID].append('{}\tpseudogene\tCDS\t{}\t{}\t.\t{}\t.\tID={}.CDS.{}\n'.format(sublist[0], str(sublist[1]), str(sublist[2]), sublist[3], dipgene, counter))
		counter += 1 

finallines = []
geneorder = list(results.keys())
with open('/global/cscratch1/sd/vstartag/RNAseqBhyb26/BhybridumBhyb26v2.1.gene.gff3', 'r') as gff:
	read = gff.readlines()
	for headerline in read[0:3]:
		finallines.append(headerline)
	currentkey = ''
	for n in range(3, len(read)):
		line = read[n]
		L = line.split('\t')
		if L[2] == 'gene':
			keyID = '{}_{}'.format(L[0], L[3])
			geneorder.append(keyID)
			results[keyID] = [line]
			currentkey = keyID
		else:
			results[currentkey].append(line)

			
geneorder = natural_sort(geneorder)
#testL = ['BhD1_900', 'BhS10_600', 'BhS1_200', 'BhD1_500', 'BhS1_250', 'BhS10_400']
for keyID in geneorder:
	for line in results[keyID]:
		finallines.append(line)

with open('/global/cscratch1/sd/vstartag/RNAseqBhyb26/Bhyb26gff_w_pseudogenes.gff3', 'w') as outF:
	for line in finallines:
		outF.write(line)

