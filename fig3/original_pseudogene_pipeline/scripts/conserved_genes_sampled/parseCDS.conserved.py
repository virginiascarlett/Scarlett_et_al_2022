import sys

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

mydir = sys.argv[1]

superCDSdict = {}
CDSdicts = []
cdsfilenames = ['Bdistachyon_556_v3.2.cds_primaryTranscriptOnly.fa', 'Bstacei_316_v1.1.cds_primaryTranscriptOnly.fa']
for f in cdsfilenames:
	CDSdicts.append(read_fasta(f))
for d in CDSdicts:
	for k, v in d.items():
		superCDSdict[k] = v

mygenes = []
with open('{}/length_differences.txt'.format(mydir), 'r') as myfile:
	for line in myfile.readlines():
		mygenes.append(line.split('\t')[0])

aliases = {}
with open('{}/slurm_aliases.txt'.format(mydir), 'r') as myfile:
	for line in myfile.readlines():
		alias, geneID = line.split('\t')[0], line.split('\t')[1].strip('\n')
		aliases[geneID] = alias

for geneID in mygenes:
	with open('{}/candidatePseudogeneRegions/{}.cds.fa'.format(mydir, aliases[geneID]), 'w') as outF:
		outF.write('>{}\n'.format(geneID))
		outF.write('{}\n'.format(superCDSdict[geneID]))


