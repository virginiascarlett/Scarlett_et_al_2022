import glob
import timeForTE

files = glob.glob('LTRs_by_family/*.trimmed2')
for f in files:
	filedic = timeForTE.read_fasta(f)
	outF = open(f, 'w')
	for header, seq in filedic.items():
		outF.write('>'+header.strip('>_R_')+'\n')
		outF.write('{}\n'.format(seq))

