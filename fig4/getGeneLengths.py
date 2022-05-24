import glob

def uniquify(seq):
	seen = set()
	seen_add = seen.add
	return [x for x in seq if not (x in seen or seen_add(x))]

def read_fasta(filename): 
	lines = {}
	with open(filename, 'r') as fastafile:
		currentKey = ''
		for line in fastafile.readlines():
			if line.startswith('>'):
				currentKey = line.strip('\n')#.strip('>')
			else:
				if currentKey in lines:
					lines[currentKey] += line.strip('\n')
				else:
					lines[currentKey] = line.strip('\n')
	return(lines)

allgenes = read_fasta('BhybridumBhyb26v2.1.primaryTrs.fa')
#For old ABR113 data:
#allgenes = read_fasta('/global/cscratch1/sd/vstartag/RNA_Seq_2018/refGenomes/Bhybridum_463_v1.1.transcript_primaryTranscriptOnly.fa')
genelens = {}
for k, v in allgenes.items():
	ID = k.split('locus=')[1].split()[0]
	seqlen = len(v)
	genelens[ID] = seqlen

gene2subg = {}
genesordered = []
#For old ABR113 data:
#with open('/global/cscratch1/sd/vstartag/RNA_Seq_2018/refGenomes/Bhybridum_plus_plastid.gff3', 'r') as myfile:
with open('BhybridumBhyb26v2.1.gene_exons.gff3', 'r') as myfile:
	read = myfile.readlines()
	for n in range(6, len(read)):
		line = read[n]
		if line.split('\t')[2] == 'gene' and not line.split('\t')[0].startswith('LT'):
			geneID = line.split('\t')[8].split('Name=')[1].strip('\n')
			subg = line.split('\t')[0][:3]
			gene2subg[geneID] = subg
			genesordered.append(geneID)

genesordered = uniquify(genesordered)

with open('BhDprimTrlengths.txt', 'w') as outF:
	for gene in genesordered:
		if gene2subg[gene] == 'BhD':
			outF.write('{}\t{}\n'.format(gene, genelens[gene]))

with open('BhSprimTrlengths.txt', 'w') as outF:
	for gene in genesordered:
		if gene2subg[gene] == 'BhS':
			outF.write('{}\t{}\n'.format(gene, genelens[gene]))		

	