"""
Parse a big file of genomic regions so that each genomic region gets its own fasta file.
"""
import os
import shutil
import glob

if os.path.exists('candidatePseudogeneRegions'):  
    shutil.rmtree('candidatePseudogeneRegions')
os.mkdir('candidatePseudogeneRegions')

fafilelist = glob.glob('*_gene_intervals.fa') #can be missing_gene_intervals.fa or conserved_gene_intervals.fa
#In missing_gene_intervals.fa, the header is the diploid gene ID
#In conserved_gene_intervals.fa, the header is the polyploid geneID
#I did this because we need the IDs to be unique

if len(fafilelist) > 1:
    raise Exception("You must have only one file ending with _gene_intervals.fa in your working directory.")

fafilename = fafilelist[0]
fafileobj = open(fafilename, 'r')
read = fafileobj.readlines()
fafileobj.close()

for n in range(len(read)):
    line = read[n]
    if line.startswith('>'):
        geneID = line[1:].strip('\n')
        outF = open('candidatePseudogeneRegions/{}.region.fa'.format(geneID), 'w')
        if fafilename == 'missing_gene_intervals.fa':
            outF.write('>PPLocSyntenicTo_{}\n'.format(geneID)) #PP stands for polyploid
        elif fafilename == 'conserved_gene_intervals.fa':
            outF.write('>{}\n'.format(geneID))
        outF.write(read[n+1])
        outF.close()

