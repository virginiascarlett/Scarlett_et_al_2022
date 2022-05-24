import os
import glob
import re
import sys

#THIS SCRIPT ONLY WORKS IF YOU HAVE NO EXTRA CUSTOM OUTPUT
#SO NO SUGAR, CIGAR, OR VULGAR, AND NO --ryo OR --showquerygff OR --showtargetgff
#UPDATE^ changed so that this accommodates ONLY my particular --ryo format

#run on command line like so:
#for fname in $( ls candidatePseudogeneRegions/*.aln.p2g ); do python3 PTCs_from_exonerate.py $fname ; done

# def natural_sort(l): 
#     convert = lambda text: int(text) if text.isdigit() else text.lower()
#     alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
#     return sorted(l, key=alphanum_key)


fname = sys.argv[1] #looks like e.g. candidatePseudogeneRegions/1.aln.p2g OR sample1/candidatePseudogeneRegions/1.aln.p2g OR candidatePseudogeneRegions/1_macse_NT.fasta
slurmID = int(''.join(filter(str.isdigit, fname.split('/')[-1].split('.')[0])))


result = []
with open(sys.argv[1], 'r') as myfile:
	read = myfile.readlines()
	if len(read) == 3: #no alignment could be made
		sys.exit()
	else:
		#remove some of the junk
		lines = [ line for line in read if line != '\n' ]
		lines = [ line for line in lines if line != '-- completed exonerate analysis\n' ]
		#NEW: remove fasta stuff at end of file
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


#destdir='/'.join(fname.split('/')[:-1]) 
destdir=''
if len(fname.split('/')) == 3:  #if this is a run of control genes
    destdir = '{}/PTCs_from_exonerate.txt'.format(fname.split('/')[0])
else:
    destdir = 'PTCs_from_exonerate.txt'

if result:
    if os.path.exists(destdir): 
        with open(destdir, 'a') as PTCfile:
            PTCfile.write('{}\t{}\t{}\n'.format(slurmID, result[0], result[1]) )
    else:
        with open(destdir, 'w') as PTCfile:
            PTCfile.write('alias\tPTC_in_seq1\tPTC_in_seq2\n')
            PTCfile.write('{}\t{}\t{}\n'.format(slurmID, result[0], result[1]) )
