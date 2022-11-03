"""
A script to choose a sampling of orthogroups from all the orthogroups
where all genomes were represented. Chooses n genes from BhD and n' genes from BhS
where those are the numbers of alignments we got from trying to pick up 
fragments of 'lost genes' using exonerate. This script outputs bed files for 
those Bhyb26 'regions' (actually they are just the coordinates of the appropriate 
Bhyb26 gene, since in the control case we actually have genes already, but this
allows my workflow to stay similar to what I did for the 'lost genes'). Also outputs
peptide files for the primary peptide of the diploid gene.

Run with three arguments: First is the analysis directory, e.g. sample1 
Next, how many BhD control genes and how many BhS control genes
will make up the sample pool for each trial? Do BhD first, BhS second. e.g.

python3 choose_conserved_controlsNEW.py sample1 224 240

(See notebook entry for Dec. 27 2021 for a little more info)

creates e.g. 
sample1/conserved_gene_intervals.bed
sample1/slurm_aliases.txt
sample1/candidatePseudogeneRegions/*.pep.fa (many)
"""
import random
import sys
import os
import shutil

mydir, numDgenes, numSgenes = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])

if os.path.exists(mydir):  
	shutil.rmtree(mydir)

os.mkdir(mydir)
os.mkdir('{}/candidatePseudogeneRegions'.format(mydir))


def read_fasta(filename): 
	lines = {}
	with open(filename, 'r') as fastafile:
		currentKey = ''
		for line in fastafile.readlines():
			if line.startswith('>'):
				currentKey = line.split()[3][6:] #get locus ID
			else:
				if currentKey in lines:
					lines[currentKey] += line.strip('\n')
				else:
					lines[currentKey] = line.strip('\n')
	return(lines)


Bd = [] #diploid gene IDs for this pool
BhD = [] #orthologous Bhyb26 geneIDs, in order
Bs = []
BhS = []
counter = 0
with open('conservedgenes.txt', 'r') as inF: #see my script bin_orthogroups.R
	read = inF.readlines()
	choices = random.sample(range(2, len(read)), numDgenes+numSgenes) #draw e.g. 464 'control genes' because we have 464 'lost genes'
	for n in choices:
		L = read[n].strip('\n').split('\t')
		if counter < numDgenes: #e.g. 224 genes from BhD and 240 genes from BhS
			Bd.append(L[1])
			BhD.append(L[2])
			counter += 1
		else:
			Bs.append(L[4])
			BhS.append(L[3])


#get a genomic region for Bhyb26, get peptide for diploids

Bhyb26gffdata = {}
with open('BhybridumBhyb26v2.1.gene.gff3', 'r') as Bhyb26gff:
	read = Bhyb26gff.readlines()
	for n in range(3, len(read)):
		L = read[n].strip('\n').split('\t')
		if L[2] == 'gene':
			geneID = L[8].split('Name=')[1].split(';')[0].strip('\n')
			Bhyb26gffdata[geneID] = [L[0], L[3], L[4]] #chrom, start, end

# with open('{}/conserved_gene_intervals.BhD.bed'.format(mydir), 'w') as outF:
# 	for n in range(len(Bd)):
# 		chrom, start, end = Bhyb26gffdata[BhD[n]][0], Bhyb26gffdata[BhD[n]][1], Bhyb26gffdata[BhD[n]][2]
# 		outF.write('{}\t{}\t{}\t{}\n'.format(chrom, start, end, Bd[n]))

# with open('{}/conserved_gene_intervals.BhS.bed'.format(mydir), 'w') as outF:
# 	for n in range(len(Bs)):
# 		chrom, start, end = Bhyb26gffdata[BhS[n]][0], Bhyb26gffdata[BhS[n]][1], Bhyb26gffdata[BhS[n]][2]
# 		outF.write('{}\t{}\t{}\t{}\n'.format(chrom, start, end, Bs[n]))	

final_dipgenes = Bd+Bs
final_Bhyb26genes = BhD+BhS

with open('{}/conserved_gene_intervals.bed'.format(mydir), 'w') as outF:
	for n in range(len(final_dipgenes)):
		chrom, start, end = Bhyb26gffdata[final_Bhyb26genes[n]][0], Bhyb26gffdata[final_Bhyb26genes[n]][1], Bhyb26gffdata[final_Bhyb26genes[n]][2]
		outF.write('{}\t{}\t{}\t{}\n'.format(chrom, start, end, final_dipgenes[n]))


Bdpeps = read_fasta('Bdistachyon_556_v3.2.protein_primaryTranscriptOnly.fa')
Bspeps = read_fasta('Bstacei_316_v1.1.protein_primaryTranscriptOnly.fa')
allpeps = Bdpeps
allpeps.update(Bspeps)

#### NEW
with open('{}/slurm_aliases.txt'.format(mydir), 'w') as slurmout:
	for n in range(len(final_dipgenes)):
		slurmout.write('{}\t{}\t{}\n'.format(n+1, final_dipgenes[n], final_Bhyb26genes[n]))

### Modified
for n in range(len(final_dipgenes)):
	with open('{}/candidatePseudogeneRegions/{}.pep.fa'.format(mydir, n), 'w') as outF:
		outF.write('>{}\n'.format(final_dipgenes[n]))
		outF.write(allpeps[final_dipgenes[n]])

