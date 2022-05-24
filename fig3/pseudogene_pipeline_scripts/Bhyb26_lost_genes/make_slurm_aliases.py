#run this script from dir ABOVE candidatePseudogeneRegions

import glob
import os

fafiles = glob.glob('candidatePseudogeneRegions/*.region.fa')

with open('slurm_aliases.txt', 'w') as outF:
	for n in range(len(fafiles)):
		geneID = fafiles[n][27:-10]
		outF.write('{}\t{}\n'.format(str(n+1), geneID))
		os.rename('candidatePseudogeneRegions/{}.region.fa'.format(geneID), 'candidatePseudogeneRegions/{}.region.fa'.format(str(n+1)))
		os.rename('candidatePseudogeneRegions/{}.pep.fa'.format(geneID), 'candidatePseudogeneRegions/{}.pep.fa'.format(str(n+1)))
