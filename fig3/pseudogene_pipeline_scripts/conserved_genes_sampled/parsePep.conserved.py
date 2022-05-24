import sys

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

mydir = sys.argv[1]

Bdpeps = read_fasta('Bdistachyon_556_v3.2.protein_primaryTranscriptOnly.fa')
Bspeps = read_fasta('Bstacei_316_v1.1.protein_primaryTranscriptOnly.fa')
allpeps = Bdpeps
allpeps.update(Bspeps)

GOIs = []
with open('{}/conserved_gene_intervals.sorted.bed'.format(mydir), 'r') as bedfile:
	for line in bedfile.readlines():
		GOIs.append(line.split("\t")[3].strip("\n"))

slurmIDs = {} #keys are geneIDs, values are 'aliases' (a number)
with open('{}/slurm_aliases.txt'.format(mydir), 'r') as slurmfile:
    for line in slurmfile.readlines():
        slurmID, gene1, gene2 = line.split('\t')[0], line.split('\t')[1], line.split('\t')[2].strip('\n')
        slurmIDs[gene1] = slurmID

for GOI in GOIs:
	outF = open('{}/candidatePseudogeneRegions/{}.pep.fa'.format(mydir, slurmIDs[GOI]), 'w')
	outF.write('>{}\n'.format(GOI))
	outF.write(allpeps[GOI])

