import glob
import sys

mydir = sys.argv[1]
fileList = glob.glob('{}/candidatePseudogeneRegions/*readyForPAML.fasta'.format(mydir))
myfile = fileList[0]
f = open(myfile, 'r')
read = f.readlines()
f.close()
geneID = myfile.split('.')[0]
ctlFile = open('{}/yn00.ctl'.format(mydir), 'w')
ctlFile.write('seqfile = {}\n'.format(myfile))
ctlFile.write('outfile = {}.yn\n'.format(geneID))
#ctlFile.write('runmode = -2\n') these are for codeml, not yn00
#ctlFile.write('CodonFreq = 2\n')
ctlFile.write('verbose = 0\n')
ctlFile.write('icode = 0\n')
ctlFile.write('weighting = 0\n')
ctlFile.write('commonf3x4 = 0\n')
ctlFile.close()



