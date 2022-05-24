"""
A script to remove LTRs that don't align well to
the rest of the family. These are poor RepeatMasker
hits, and may or may not even be real LTRs, and they will
screw up our trees. 

We have already obtained MSAs for the LTR families using
MAFFT and TrimAl. Now I am using those MSAs to determine 
sequence quality. I am removing LTR sequences that 
are more than 20% gaps, creating a new fasta file called
*.cleaned.fa. This file is identical to the original
fasta file that was used as MAFFT input, except that 
questionable fasta entries have been removed. This cleaned 
file will be the input for a new round of MAFFT + TrimAl,
a do-over basically.

I am creating the batch script that will run a job 
array of MAFFT + TrimAl jobs. I am also creating the batch 
scripts required to run raxml, which we will run once 
mafft + trimal is finished. Note that some LTR families have 
a .constraint file, which I am creating here.
This will tell raxml that 5' and 3' mates of a single fl-LTR-RT
should be in the same monophyletic group (though there might be 
additional sequences in the group, too). Families with no fl-LTR-RTs, 
i.e. comprised of only solo-LTRs, will not have a constraint 
file, since we have no a priori knowledge about their topology.

Note that we will have cases where one LTR from a fl-LTR-RT
is discarded but the other is retained. The remaining LTR will 
be treated as a solo LTR for the purposes of building constraint trees.
In the next script, raxml_terminal_branch_lengths.py, we will 
label this TE as not an fl-LTR-RT. 

EDIT LINES 194 AND 211 IF YOU WANT TO USE A MODEL OTHER THAN GTR+G4!!!
"""

import glob
import more_itertools as mit
import timeForTE


def print_slurm_codes(finalDict):
	slurmcodes = sorted([int(x) for x in finalDict.keys()])
	slurm_codes_grouped = [list(group) for group in mit.consecutive_groups(slurmcodes)]
	finalcodes = []
	for sublist in slurm_codes_grouped:
		if len(sublist)==1:
			finalcodes.append(str(sublist[0]))
		elif len(sublist)>1:
			finalcodes.append('-'.join((str(sublist[0]), str(sublist[-1]))))
	return( ','.join(finalcodes) )

def check_for_discarded_mates(flLTRRTs_frags):
	onlyOneMate = []
	for TEcode, fraglist in flLTRRTs_frags.items():
		if len(fraglist) != 2:
			onlyOneMate.append(TEcode)
			if len(fraglist) != 1: #bugtesting
				print('ERROR! ERROR! At least one fl-LTR-RT in the family with slurmID {} has either 0 or >2 fragments!'.format(slurmID))
	return(onlyOneMate)

def get_trimmed_dic(filename):
	trimmed_seqs = {}
	read = open(filename, 'r').readlines()
	frags = []
	data = []
	current_frag_data = []
	for line in read:
		if line.startswith('>'):
			if frags:
				data.append(''.join(current_frag_data))
			current_frag_data = []
			frags.append('>'+line.strip('>_R_'))
		elif not line.startswith('>'):
			current_frag_data.append(line.rstrip())
	data.append(''.join(current_frag_data))
	for i in range(len(frags)):
		trimmed_seqs[frags[i]] = data[i]
	return(trimmed_seqs)


def get_raw_dic(filename):
	raw_seqs = {}
	originalfasta = open(filename, 'r').readlines()
	for n in range(len(originalfasta)):
		if originalfasta[n].startswith('>'):
			raw_seqs[originalfasta[n]] = originalfasta[n+1]
	return(raw_seqs)

def get_keepers(trimdic): #returns a list of fragments worth keeping
	to_keep = []
	for header, seq in trimdic.items():
		num_gaps = trimdic[header].count('-')
		seqlen = len(trimdic[header])
		if num_gaps/seqlen < threshold: 
			to_keep.append(header)
	if len(to_keep) < 4: #Need at least 4 sequences to make a tree
		to_keep = []
	return(to_keep)





files = glob.glob('LTRs_by_family/*.trimmed')
threshold = 0.20

seqsByFamily = {} #a dict with key=family name and value=a dict.
#In the subdict, key=a fragment code (or fragcodes concatenated with '_'), value=the DNA sequence
for f in files:
	slurmID = f.split('/')[1].replace('.trimmed', '')
	seqsByFamily[slurmID] = {}
	trimmed_seqs = get_trimmed_dic(f)
	seqsByFamily[slurmID] = trimmed_seqs


originalFastaData = {}
for slurmID in seqsByFamily.keys():
	rawfilename = 'LTRs_by_family/{}.fa'.format(slurmID)
	raw_seqs = get_raw_dic(rawfilename)
	originalfasta = open('LTRs_by_family/{}.fa'.format(slurmID), 'r').readlines()
	originalFastaData[slurmID] = raw_seqs


seqsByFamilyFinal = {}  #e.g. list(seqsByFamilyFinal['57'].keys()) yields
#['>frag0006544\n', '>frag0006546\n', '>frag0016437\n', '>frag0016439\n', '>frag0055159\n', '>frag0056790\n', '>frag0056795\n', '>frag0061462_frag0061463\n']
for slurmID, fastadict in seqsByFamily.items():
	keepers = get_keepers(fastadict)
	if keepers:
		seqsByFamilyFinal[slurmID] = {}
		for header in keepers:
			seqsByFamilyFinal[slurmID][header] = originalFastaData[slurmID][header]


for slurmID, outdic in seqsByFamilyFinal.items():
	with open('LTRs_by_family/{}.cleaned.fa'.format(slurmID), 'w') as outF:
		for header, seq in outdic.items():
			outF.write(header)
			outF.write(seq)

#Create constraint files where applicable
hella_data_structures = timeForTE.read_fragments_classified()
fragstocopies = hella_data_structures[3] #key is fragcode and value is TEcopycode
flLTRRTs = hella_data_structures[5] #indexed by TEcopy code, value is TEcopy object

seqsByFamilyConstraint = {} #families with at least one fl-LTR-RT
seqsByFamilyNoConstraint = {} #families with no fl-LTR-RTs
quest = open('LTRs_by_family/questionableLTRs.txt', 'w')
quest.write('flLTRRT\tremainingLTR\n')
for slurmID, outdic in seqsByFamilyFinal.items():
	flLTRRTs_frags = {} #frags in this family that belong to an fl-LTR-RT. key=TEcode, value=list containing two fragcodes
	remaining_frags = []
	for header in outdic.keys():
		myfragcode = header.replace('>', '').replace('\n', '')
		myTEcode = fragstocopies[myfragcode.split('_')[0]]
		if myTEcode in flLTRRTs:
			if myTEcode in flLTRRTs_frags:
				flLTRRTs_frags[myTEcode].append(myfragcode)
			else:
				flLTRRTs_frags[myTEcode] = [myfragcode]
		else:
			remaining_frags.append(myfragcode)
	if flLTRRTs_frags:
		toRemove = check_for_discarded_mates(flLTRRTs_frags)
		if toRemove:
			for TEcode in toRemove:
				remaining_frags.append(flLTRRTs_frags[TEcode][0])
				quest.write('{}\t{}\n'.format(TEcode, flLTRRTs_frags[TEcode][0]))
				del(flLTRRTs_frags[TEcode])
		seqsByFamilyConstraint[slurmID] = outdic
		flLTRstr = ','.join([ str((fragpair[0], fragpair[1])) for fragpair in flLTRRTs_frags.values()]).replace(' ', '').replace("'", "")
		remainingstr = str( ','.join([h.replace('>', '').replace('\n', '') for h in  remaining_frags]) )
		constraintFileName = 'LTRs_by_family/{}.constraint'.format(slurmID)
		with open(constraintFileName, 'w') as conOut:
			conOut.write('({},{})\n'.format(remainingstr, flLTRstr))
	else:
		seqsByFamilyNoConstraint[slurmID] = outdic



#write_SBATCH(nodes, ntasks, cpus_per_task, time, mem, qos, mailtype, filename)
timeForTE.write_SBATCH('1', '1', '32', '00:45:00', '59G', timeForTE.qos_partial_node, 'fail', 'LTRs_by_family/runmafft2.sh')
with open('LTRs_by_family/runmafft2.sh', 'a') as batcharray:
	batcharray.write('#SBATCH --array={}\n'.format(print_slurm_codes(seqsByFamilyFinal)))
	batcharray.write('#SBATCH --error=align2.%A.%a.err\n')
	batcharray.write('#SBATCH --output=align2%A.%a.out\n')
	batcharray.write('\n')
	batcharray.write('einsi --thread -1 --adjustdirectionaccurately ${SLURM_ARRAY_TASK_ID}.cleaned.fa > ${SLURM_ARRAY_TASK_ID}.aln2\n')
	batcharray.write('sleep 2\n')
	batcharray.write('trimal -in ${SLURM_ARRAY_TASK_ID}.aln2 -gt 0.8 -out ${SLURM_ARRAY_TASK_ID}.trimmed2\n')

timeForTE.write_SBATCH('1', '1', '1', '3:00:00', '2G', timeForTE.qos_partial_node, 'fail', 'LTRs_by_family/runraxml.noconstraint.sh')
with open('LTRs_by_family/runraxml.noconstraint.sh', 'a') as batcharray:
	batcharray.write('#SBATCH --array={}\n'.format(print_slurm_codes(seqsByFamilyNoConstraint)))
	batcharray.write('#SBATCH --error=tree.%A.%a.err\n')
	batcharray.write('#SBATCH --output=tree.%A.%a.out\n')
	batcharray.write('\n')
	batcharray.write('raxml-ng --msa ${SLURM_ARRAY_TASK_ID}.trimmed2 --model GTR+G4 --threads 1 --redo\n')


timeForTE.write_SBATCH('1', '1', '1', '1:00:00', '2G', timeForTE.qos_partial_node, 'fail', 'LTRs_by_family/runraxml.constraint.sh')
with open('LTRs_by_family/runraxml.constraint.sh', 'a') as batcharray:
	batcharray.write('#SBATCH --array={}\n'.format(print_slurm_codes(seqsByFamilyConstraint)))
	batcharray.write('#SBATCH --error=tree.%A.%a.err\n')
	batcharray.write('#SBATCH --output=tree.%A.%a.out\n')
	batcharray.write('\n')
	batcharray.write('raxml-ng --msa ${SLURM_ARRAY_TASK_ID}.trimmed2 --model GTR+G4 --threads 1 --redo --tree-constraint ${SLURM_ARRAY_TASK_ID}.constraint\n')


#SCRAPS 
# >>> list(trimmed_seqs.keys())[0:5]
# ['>frag0000086\n', '>frag0000089\n', '>frag0001271\n', '>frag0001273\n', '>frag0003874\n']
# >>> list(trimmed_seqs.values())[0:5]
# ['cagaaatatgagatta--------------taatttttatttatctgaccctcgttggtatatataatacgagggataaactggtacaagacaagtaatttagattctatcttatctcccata', 'tagaaatatgagattgatgcttaattgatatatgttgtattgatctgatcctcgttgatatatatgatacagagcatagactaatacaagacaattaatttagattctatc------------', 'tagagaaacaagatagaattataaacattatgtcttgtactagtatattccccctttatatatataatatgagggtcagatcag---aacacaaatactttaggttcaatcccatattctaca', 'tagagagataagatagaattataaacattatgtcttgtactagtatatcccccctttgtatatataatatgaggatcagatcaa---aacacaaatactttaggttcaatcccatattctaca', 'tggaaatatgagattgaagcctaattggtatatgttgtattgatctgaccctcattggtatatatagtacgggggatagactggtacaagacaagtaatttaaattctatcttatctcttata']


# >>> list(raw_seqs.keys())[0:5]
# ['>frag0000086\n', '>frag0000089\n', '>frag0001271\n', '>frag0001273\n', '>frag0003874\n']
# >>> list(raw_seqs.values())[0:5]
# ['GCAGAAATATGAGATTATTTCACGGTATACGGTGATTTTTATTTATCTGACCCTCGTATGGGGTATATATAGAGATACAATGGAGGAGATAAACTTGGAGTACAAGACAAGTAATTTCTTTAAACCGATTCTATCTTATCTCTAATCCCAATCATACCATAGCTCTAACA\n', 'GTAGAAATATGAGATTGGATGCTTAATTGATATTCCACGATATACTGTGTTGTTGTATTGATCTGATCCTCGTATTGGATATATATGGAGATACAATGAGAGACATAGACTTAAAATACAAGACAATTAATCTTTTTAAACCGATTCTATC\n', 'GTTAGAGATATAATACGGTTAGGACTAGAAACAAGATAGAATTAGTTTAAACGAGATTACTTGTCTTGTACTCCAAGTATATTTCCCCCTTTATACCTCTATATATACCCCATATGAGGGTCAGATCAGAACACAGTATATCGTGGAATACAATTTAGGGTTCCAATCCCATATTTCTACA\n', 'GTTAGAGATATAATACGGTTAGGATTAGAGATAAGATAGAATTAGTTTAAACGAGATTACTTGTCTTGTACTCCAAGTATATCTCCCCCTTTGTACCTCTATATATACCCCATATGAGGATCAGATCAAAACACAGTATATCGTGGAATACAATTTAGGGTTCCAATCCCATATTTCTACA\n', 'GTGGAAATATGAGATTGGAAGCCTAATTGGTATTCCACGATATTCGTTGTTGTATTGATCTGACCCTCATATGGGGTATATATAGAGGTACAAAGGGGGAGATAGACTTGGAGTACAAGACAAGTAATCTTGTTTAAACTAATTCTATCTTATCTCTAATCCCAAACATATTATATCTCTAACA\n']


# >>> to_keep[0:5]
# ['>frag0000086\n', '>frag0000089\n', '>frag0001271\n', '>frag0001273\n', '>frag0003874\n']


