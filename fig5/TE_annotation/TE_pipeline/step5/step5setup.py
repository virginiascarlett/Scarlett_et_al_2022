"""
A script to take the per-chromosome output files produced by 
RepeatMasker in step4 and concatenate them into a single gff.
Run from your step5/ dir. The files should be in your step4/ dir.
"""

import timeForTE
import glob

alloutfiles = glob.glob('{}/step4/*.out'.format(timeForTE.topdir))
outfiles = []
#make sure these are actually RepeatMasker output files
for fname in alloutfiles:
	with open(fname, 'r') as myfile:
		if myfile.readlines()[0].strip().startswith('SW'):
			outfiles.append(fname)


outfiles.sort(key=lambda name: int(name.split('/')[-1].split('.')[0]))
outlines = []
for fname in outfiles:
	with open(fname, 'r') as fobj:
		read = fobj.readlines()
		if outlines == []: #if this is the first file we're processing, keep the header
			outlines = [ line for line in read ]
		else:
			for n in range(3, len(read)): #skip header
				outlines.append(read[n])


outF = open( '{}.all.fasta.out'.format(timeForTE.genome), 'w' )			
for line in outlines:
	outF.write(line)

outF.close()
