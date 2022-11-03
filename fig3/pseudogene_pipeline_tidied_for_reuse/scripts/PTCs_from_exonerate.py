"""
A script to count how many premature termination codons (PTCs) occur in an
exonerate alignment. Produces an output file called PTCs_from_exonerate.txt.
THIS SCRIPT ONLY WORKS IF YOU RAN EXONERATE WITH A PARTICULAR --ryo format:
--ryo ">%qi %ql %qal\n%qas\n>%ti %tcl\n%tcs\n"

run this script like so:
for fname in $( ls candidatePseudogeneRegions/*.aln.p2g ); do python3 PTCs_from_exonerate.py $fname ; done

In the output file, 'seq1' refers to the first sequence (query) in the exonerate alignment
and 'seq2' refers to the second sequence (target).
"""


import os
import glob
import re
import sys

fname = sys.argv[1] #looks like e.g. candidatePseudogeneRegions/1.aln.p2g 
alias = int(''.join(filter(str.isdigit, fname.split('/')[-1].split('.')[0])))


result = []
with open(fname, 'r') as myfile:
	read = myfile.readlines()
	if len(read) == 3: #no alignment could be made
		sys.exit()
	else:
		#remove some of the junk
		lines = [ line for line in read if line != '\n' ]
		lines = [ line for line in lines if line != '-- completed exonerate analysis\n' ]
		#remove fasta stuff at end of file
		myaln = []
		for n in range(len(lines)):
			if lines[n].startswith('>'):
				break
			else:
				myaln.append(lines[n])
		pattern = re.compile(r'\s+') #matches all whitespace (space, tab, newline, and so on). we'll remove this.
		queryAA = ''
		for i in range(10, len(myaln), 4):
			queryAA += re.sub(pattern, '', myaln[i]).split(':')[1]
		targetAA = ''
		for i in range(12, len(myaln), 4):
			targetAA += re.sub(pattern, '', myaln[i])
		#remove the real termination codon if there is one
		if queryAA.endswith('***') and targetAA.endswith('***'):
			queryAA = queryAA[:-3]
			targetAA = targetAA[:-3]
		# I'm just recording whether '***' occurs, (True/False) not the number of occurrences.
		if '***' in queryAA or '***' in targetAA:
			result = ['***' in queryAA, '***' in targetAA]

if result:
	if os.path.exists('PTCs_from_exonerate.txt'): 
		with open('PTCs_from_exonerate.txt', 'a') as PTCfile:
			PTCfile.write('{}\t{}\t{}\n'.format(alias, result[0], result[1]) )
	else:
		with open('PTCs_from_exonerate.txt', 'w') as PTCfile:
			PTCfile.write('alias\tPTC_in_seq1\tPTC_in_seq2\n')
			PTCfile.write('{}\t{}\t{}\n'.format(alias, result[0], result[1]) )

