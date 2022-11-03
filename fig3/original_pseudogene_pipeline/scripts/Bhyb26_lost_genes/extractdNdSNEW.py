"""
For control genes, run like this:
python3 extractdNdS.conserved.py sample${SLURM_ARRAY_TASK_ID}
For lost genes, don't include any command line args
"""
import glob
import sys

mydir = ''
if len(sys.argv) > 1:
    mydir = sys.argv[1]+'/'

slurmIDs2genes = {}
with open('{}slurm_aliases.txt'.format(mydir), 'r') as slurmfile:
    for line in slurmfile.readlines():
        alias, gene = line.split('\t')[0], line.split('\t')[1].strip('\n')
        slurmIDs2genes[alias] = gene


fileList = glob.glob('{}candidatePseudogeneRegions/*.yn'.format(mydir))
outF = open('{}dNdSbyGene.txt'.format(mydir), 'w')
bugtest = []
for fname in fileList:
    f = open(fname, 'r')
    read = f.readlines()
    f.close()
    if len(read) > 0:
        geneID = slurmIDs2genes[ fname.strip('.yn')[27:] ]
        info = list(filter(lambda a: a != '', read[90].split(' '))) #parse the 90th line, which has the stuff we want
        omega = info[6]
        outF.write('{}\t{}\n'.format(geneID, omega))
    else:
        bugtest.append(fname)

outF.close()
