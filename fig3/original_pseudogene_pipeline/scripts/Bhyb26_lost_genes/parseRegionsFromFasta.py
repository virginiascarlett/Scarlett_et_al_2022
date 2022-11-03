"""
Parse a big file of genomic regions created by identifyTriplets.py followed by bedtools getfasta.
Split up so that each genomic region gets its own fasta file.
"""
import os
import shutil

if os.path.exists('candidatePseudogeneRegions'):  
    shutil.rmtree('candidatePseudogeneRegions')
os.mkdir('candidatePseudogeneRegions')

fafile = open('missing_gene_intervals.fa', 'r')

read = fafile.readlines()
for n in range(len(read)):
    line = read[n]
    if line.startswith('>'):
        geneID = line[1:].strip('\n')
        outF = open('candidatePseudogeneRegions/{}.region.fa'.format(geneID), 'w')
        outF.write('>BhybLocSyntenicTo_{}\n'.format(geneID))
        outF.write(read[n+1])
        outF.close()

fafile.close()


