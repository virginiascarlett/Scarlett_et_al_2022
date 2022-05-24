"""
A script to get Bhyb26 genomic coordinates for all widely conserved Bhyb26 genes. 
We will ultimately get dN/dS values for these genes relative to their diploid orthologs. 
Actually, we are just taking all fully conserved syntenic orthogroups and arbitrarily 
(50/50 chance) picking either the BhD or BhS gene (15,000 genes sould be enough--30,000 
is just a lot more I/O and seems overkill).

We are also making pep.fa files for the diploid genes.

Just run once with no command line args. Note that if candidatePseudogeneRegions/
exists, this script will remove it.
"""
import random
import os
import shutil

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


if os.path.exists('candidatePseudogeneRegions'):  
	shutil.rmtree('candidatePseudogeneRegions')

os.mkdir('candidatePseudogeneRegions')

#get a genomic region for Bhyb26
Bhyb26gffdata = {}
with open('BhybridumBhyb26v2.1.gene.gff3', 'r') as Bhyb26gff:
	read = Bhyb26gff.readlines()
	for n in range(3, len(read)):
		L = read[n].strip('\n').split('\t')
		if L[2] == 'gene':
			geneID = L[8].split('Name=')[1].split(';')[0].strip('\n')
			Bhyb26gffdata[geneID] = [L[0], L[3], L[4]] #chrom, start, end


outF = open('conserved_gene_intervals.bed', 'w')

final_dipgenes = []
final_Bhyb26genes = []
with open('conservedgenes.txt', 'r') as inF: #see my script bin_orthogroups.R
	read = inF.readlines()
	for n in range(1, len(read)):
		L = read[n].strip('\n').split('\t')
		Bdgene, Bhyb26Dgene, Bhyb26Sgene, Bsgene = L[1], L[2], L[3], L[4]
		result = random.choices(['BhD', 'BhS'], k = 1)[0]
		if result == 'BhD':
			chrom, start, end = Bhyb26gffdata[Bhyb26Dgene][0], Bhyb26gffdata[Bhyb26Dgene][1], Bhyb26gffdata[Bhyb26Dgene][2]
			outF.write('{}\t{}\t{}\t{}\n'.format(chrom, start, end, Bdgene))
			final_dipgenes.append(Bdgene)
			final_Bhyb26genes.append(Bhyb26Dgene)
		elif result == 'BhS':
			chrom, start, end = Bhyb26gffdata[Bhyb26Sgene][0], Bhyb26gffdata[Bhyb26Sgene][1], Bhyb26gffdata[Bhyb26Sgene][2]
			outF.write('{}\t{}\t{}\t{}\n'.format(chrom, start, end, Bsgene))
			final_dipgenes.append(Bsgene)
			final_Bhyb26genes.append(Bhyb26Sgene)

outF.close()

#get peptide for diploids
Bdpeps = read_fasta('Bdistachyon_556_v3.2.protein_primaryTranscriptOnly.fa')
Bspeps = read_fasta('Bstacei_316_v1.1.protein_primaryTranscriptOnly.fa')
allpeps = Bdpeps
allpeps.update(Bspeps)

slurmout = open('slurm_aliases.txt', 'w') 
for n in range(len(final_dipgenes)):
	with open('candidatePseudogeneRegions/{}.pep.fa'.format(n+1), 'w') as outF: #goddamn zero-based numbering
		mydipgene, myBhgene = final_dipgenes[n], final_Bhyb26genes[n]
		outF.write('>{}\n'.format(mydipgene))
		outF.write(allpeps[mydipgene])
		slurmout.write('{}\t{}\t{}\n'.format(n+1, mydipgene, myBhgene))

slurmout.close()

