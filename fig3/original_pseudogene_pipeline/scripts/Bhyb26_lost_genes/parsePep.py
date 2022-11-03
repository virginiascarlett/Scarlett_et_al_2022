
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

Bdpeps = read_fasta('Bdistachyon_556_v3.2.protein_primaryTranscriptOnly.fa')
Bspeps = read_fasta('Bstacei_316_v1.1.protein_primaryTranscriptOnly.fa')
allpeps = Bdpeps
allpeps.update(Bspeps)

GOIs = []
with open('missing_gene_intervals.sorted.bed', 'r') as bedfile:
	for line in bedfile.readlines():
		GOIs.append(line.split("\t")[3].strip("\n"))

for GOI in GOIs:
	outF = open('candidatePseudogeneRegions/{}.pep.fa'.format(GOI), 'w')
	outF.write('>{}\n'.format(GOI))
	outF.write(allpeps[GOI])

