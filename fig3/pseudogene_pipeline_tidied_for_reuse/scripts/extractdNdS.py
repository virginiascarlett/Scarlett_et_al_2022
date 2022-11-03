"""
For control genes, run like this:
python3 extractdNdS.conserved.py sample${SLURM_ARRAY_TASK_ID}
For lost genes, don't include any command line args
"""
import glob

fileList = glob.glob('candidatePseudogeneRegions/*.yn')

aliases = {} # If we are working with missing genes, the values are diploid genes.
#If we are working with conserved genes, the values are polyploid genes. Keys are their aliases.
with open('aliases.txt', 'r') as aliasfile:
    for line in aliasfile.readlines():
        alias, geneID = line.split('\t')[0], line.split('\t')[1].strip('\n')
        aliases[alias] = geneID

outF = open('dNdSbyGene.txt', 'w') 
for fname in fileList:
    f = open(fname, 'r')
    read = f.readlines()
    f.close()
    if len(read) > 0: #Some files will be empty, meaning an omega value couldn't be found
        alias = fname.strip('.yn')[27:]
        geneID = aliases[alias]
        info = list(filter(lambda a: a != '', read[90].split(' '))) #parse the 90th line, which has the stuff we want
        omega = info[6]
        outF.write('{}\t{}\n'.format(geneID, omega))

outF.close()
