"""
Parse a big file of genomic regions created by identifyTriplets.py followed by bedtools getfasta.
Split up so that each genomic region gets its own fasta file.
"""

slurmIDs = {} #keys are geneIDs, values are 'aliases' (a number)
gene1_to_gene2 = {} #keys are diploid genes, values are hybridum genes (or gene1 and gene2)
with open('slurm_aliases.txt', 'r') as slurmfile:
    for line in slurmfile.readlines():
        slurmID, gene1, gene2 = line.split('\t')[0], line.split('\t')[1], line.split('\t')[2].strip('\n')
        slurmIDs[gene1] = slurmID
        gene1_to_gene2[gene1] = gene2

fafile = open('conserved_gene_intervals.fa', 'r')
read = fafile.readlines()
fafile.close()
for n in range(len(read)):
    line = read[n]
    if line.startswith('>'):
        geneID = line[1:].strip('\n')
        outF = open('candidatePseudogeneRegions/{}.region.fa'.format(slurmIDs[geneID]), 'w')
        outF.write( '>{}\n'.format(gene1_to_gene2[geneID]) )
        outF.write(read[n+1])
        outF.close()


