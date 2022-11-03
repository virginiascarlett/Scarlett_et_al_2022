"""
PAML uses a config file called the 'control file'. 
This script makes one control file for the first gene/pseudogene pair.
Keeping with the PAML convention, it's called yn00.ctl.
My script run_yn00.sh iteratively updates the control file
and then runs yn00 on a gene pair, in a loop for all pairs. 
"""
import glob

fileList = glob.glob('candidatePseudogeneRegions/*.readyForPAML.fasta')
myfile = fileList[0]
f = open(myfile, 'r')
read = f.readlines()
f.close()
geneID = myfile.split('.')[0][27:]
ctlFile = open('yn00.ctl', 'w')
ctlFile.write('seqfile = {}\n'.format(myfile))
ctlFile.write('outfile = candidatePseudogeneRegions/{}.yn\n'.format(geneID))
#ctlFile.write('runmode = -2\n') these are for codeml, not yn00
#ctlFile.write('CodonFreq = 2\n')
ctlFile.write('verbose = 0\n')
ctlFile.write('icode = 0\n')
ctlFile.write('weighting = 0\n')
ctlFile.write('commonf3x4 = 0\n')
ctlFile.close()