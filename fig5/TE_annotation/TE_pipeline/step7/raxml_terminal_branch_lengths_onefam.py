"""
A script to collect the terminal branch lengths from all the raxml 
tree files in LTRs_by_family/, and puts them into one file.
Run while in your step7/ dir.
Output file has columns: familyname, TEcode, fragcode(s), terminal_branch_length.
If a particular TE has two LTRs represented in the tree, I am
summing their respective terminal branch lengths to get one final branch length for that TE.
Recall that terminal branch lengths will be an upper limit age estimate for solo-LTRs. 


For simplicity: if the input sequence prior to alignment consisted of multiple 
concatednated fragments, I am just writing the name of the first fragment. If two 
LTRs from a single fl-LTR-RT were summed to get one branch length, they will be 
joined by a ','.

Output file columns are RMname (family name), TE copy code, terminal branch length and insertion time branchlength/(2*mutation_rate)

"""

import os
import sys
import timeForTE
from Bio import Phylo

def group_mates(branchlist):
	outlist = []
	for tup in branchlist:
		myfragcode = fragobjs[tup[0].split('_')[0]].fragcode
		myTEcode = fragobjs[tup[0].split('_')[0]].TEcode
		added = 'no'
		for minidict in outlist:
			if fragobjs[list(minidict.keys())[0]].TEcode == myTEcode:
				minidict[myfragcode] = tup[1]
				added = 'yes'
		if added == 'no':
			outlist.append({myfragcode:tup[1]})
	return(outlist)

	

if len(sys.argv)!=2:
	print('Please include the slurm alias of your family of interest as a command line argument. E.g. python3 raxml_terminal_branch_lengths_onefam.py 62')
	sys.exit()

family_of_interest = sys.argv[1]
filename = '{}.trimmed2.raxml.bestTree'.format(family_of_interest)


hella_data_structures = timeForTE.read_fragments_classified()
TEcopyobjs = hella_data_structures[0] #indexed by TEcopy code
fraglists = hella_data_structures[6] #indexed by TEcopy code, a list of TEfrag objects
fragobjs = hella_data_structures[1] #key is fragcode and value is fragobj
fragstocopies = hella_data_structures[3] #key is fragcode and value is TEcopycode
flLTRRTs = hella_data_structures[5] #indexed by TEcopy code, value is TEcopy object
exobjs = hella_data_structures[2] #key is TEcopy code, value is a list: [exemplarobj, classification]

os.chdir ("LTRs_by_family")


tree = Phylo.read(filename, "newick")
branchlist = [(t.name, t.branch_length) for t in tree.get_terminals()] 
#branchlist looks like e.g. [('frag0055383', 0.052344), ('frag0010172', 0.057267), ('frag0085185', 0.040003)]
term_branch_lengths = group_mates(branchlist)
#term_branch_lengths is a list of dicts.
#if a given TE copy is represented by only one fragment in the tree, the dict will have only one entry.
#If a 5' and a 3' LTR is present, it will have two. dict key is a fragcode, dict value is terminal branch length.
#Looks like e.g.: [{'frag0024250': 0.051447, 'frag0024247': 0.078408}, {'frag0010989': 0.131079}, {'frag0063129': 0.01526, 'frag0063131': 0.027004}]


final_branch_lengths = {} #looks similar to term_branch_lengths, but 3' and 5' LTRs 
#from fl-LTR-RTs are represented as a single sum, and their names are joined with a ','. 
#Also final_branch_lengths is a flat dict, instead of a list of dicts. Looks like e.g.
# {'frag0024250,frag0024247': 0.129855, 'frag0010989': 0.131079, 'frag0063129,frag0063131': 0.04226399999}
for d in term_branch_lengths: 
	if len(d) == 2:
		newfragcode = ','.join(d.keys())
		newbranch = sum(d.values())
		final_branch_lengths[newfragcode] = newbranch
	else:
		final_branch_lengths[list(d.keys())[0]] = list(d.values())[0]

insertion_times = {} # a flat dict, looks like e.g.:
#{'frag0218907': 4638538.461538461, 'frag0058255': 3540307.6923076925, 'frag0212519,frag0212521': 1810461.5384615385}
mutation_rate = 1.3E-8
for fragcode, branchlength in final_branch_lengths.items():
	age = branchlength/(2*mutation_rate)
	insertion_times[fragcode] = age



os.chdir("..")
outF = ''
outFname = 'raxml_collected.out_{}_topfams'.format(timeForTE.genome)
if os.exists(outFname):
	outF = open(outFname, 'a')
else:
	outF = open(outFname, 'w')
	outF.write('exemplar\tTEcode\tfrag(s)\tbranchlen\tinsertion_time\tisFL\tclassification\n')

for fcode, branchlength in final_branch_lengths.items():
	myTEcode = fragstocopies[fcode.split(',')[0]]
	outF.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
		TEcopyobjs[myTEcode].RMname, 
		myTEcode, 
		fcode,
		str(branchlength),
		str(round(insertion_times[RMname][fcode], 2)),
		str(timeForTE.is_flLTR_RT(fraglists[myTEcode])),
		exobjs[myTEcode][1]
		)
	)


outF.close()
