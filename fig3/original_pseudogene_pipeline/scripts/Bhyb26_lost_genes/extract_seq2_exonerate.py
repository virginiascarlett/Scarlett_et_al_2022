"""
This script will pull the second sequence (target) from an exonerate output file
that was made with this format specification:
--ryo ">%qi %ql %qal\n%qas\n>%ti %tcl\n%tcs\n" 
and dump it into a fasta file.

This ryo specification looks like:
>Dipgene num_AAs_in_peptide num_AAs_in_alignment
M...*
>BhybLocSyntenicTo_Dipgene num_NTs_in_alignment
ATG...TGA

Include file name as a command line argument e.g.:
python3 extract_seq1_exonerate.py candidatePseudogeneRegions/1.aln.p2g

To run this on all alignment files in candidatePseudogeneRegions, do:
for fname in $( ls candidatePseudogeneRegions/*.aln.p2g ); do python3 extract_seq2_exonerate.py $fname ; done
"""
import os
import sys


class ExonerateAlignment(object):
	def __init__(self, gene1, seq1, gene2, seq2, seq1_initlen, seq1_finallen, seq2_finallen):
		self.gene1 = gene1
		self.seq1 = seq1
		self.gene2 = gene2
		self.seq2 = seq2
		self.seq1_initlen = int(seq1_initlen)
		self.seq1_finallen = int(seq1_finallen)
		self.seq2_finallen = int(seq2_finallen)


def parse_exonerate(inFname):
	with open(inFname, 'r') as myfile:
		g1id, s1, g2id, s2, s1_initlen, s1_finallen, s2_finallen = '','','','','','',''
		read = myfile.readlines()
		if len(read) == 3: #no alignment could be made
			sys.exit()
		past_seq1 = False
		for n in range(1, len(read)):
			line = read[n]
			if line.startswith(">"):
				if not past_seq1:
					g1, s1_initlen, s1_finallen = line.split()[0].strip('>'), line.split()[1], line.strip('\n').split()[2]
					counter = 1
					while read[n+counter] != '\n':
						counter +=1
					for i in range(n+1, n+counter):
						s1 += read[i].strip('\n')
					past_seq1 = True
				elif past_seq1:
					g2, s2_finallen = line.split()[0].strip('>'), line.strip('\n').split()[1]
					counter = 1
					while read[n+counter] != '\n':
						counter +=1
					for i in range(n+1, n+counter):
						s2 += read[i].strip('\n')
		myexobj = ExonerateAlignment(g1, s1, g2, s2, s1_initlen, s1_finallen, s2_finallen)
		return(myexobj)




fname = sys.argv[1] #looks like e.g. candidatePseudogeneRegions/1.aln.p2g
#OR for conserved genes looks like sample1/candidatePseudogeneRegions/1.aln.p2g
slurmID = fname.split('/')[-1][:-8]
myaln = parse_exonerate(fname)
length_diff = myaln.seq1_initlen*3 - myaln.seq1_finallen*3 

outFprefix = '/'.join(fname.split('/')[:-1])

with open('{}/{}.aln.seq2out'.format(outFprefix, slurmID), 'w') as outF:
	outF.write('>{}\n'.format(myaln.gene2))
	outF.write('{}\n'.format(myaln.seq2))


outFname = ''
if len(fname.split('/')) == 3:  #if this is a run of control genes
	outFname = '{}/length_differences.txt'.format(fname.split('/')[0])
else:
	outFname = 'length_differences.txt'


if os.path.exists(outFname):
	with open(outFname, 'a') as lenfile:
		lenfile.write('{}\t{}\n'.format(myaln.gene1, length_diff))
else:
	with open(outFname, 'w') as lenfile:
		lenfile.write('{}\t{}\n'.format(myaln.gene1, length_diff))

