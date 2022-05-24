"""
A script that takes allfamilies.classified as input and produces fasta files
for aligning all LTRs from all TE copies in a given TE family, using mafft. 
If multiple LTR sequences are present for a single TE copy, not separated by 
an int sequence, then they will be concatenated.

So for an LTR-RT that looks like this:
frag0001292     Bd1     5308156 5308404 LTR     onecode001168   5308156 5310572 +       Copia-51_BD#LTR/Copia   True    LTR_retrotransposon_RLC
frag0001293     Bd1     5308392 5308432 LTR     onecode001168   5308156 5310572 +       Copia-51_BD#LTR/Copia   True    LTR_retrotransposon_RLC
frag0001294     Bd1     5308402 5310294 int     onecode001168   5308156 5310572 +       Copia-51_BD#LTR/Copia   True    LTR_retrotransposon_RLC
frag0001295     Bd1     5310280 5310543 LTR     onecode001168   5308156 5310572 +       Copia-51_BD#LTR/Copia   True    LTR_retrotransposon_RLC
frag0001296     Bd1     5310531 5310572 LTR     onecode001168   5308156 5310572 +       Copia-51_BD#LTR/Copia   True    LTR_retrotransposon_RLC

The fasta file will contain:
>frag0001292_frag0001293
GTTGAGAAATAGCAATGCATACTTGAGAATATTCTGAAACACTCTAGATAATTGGAGAGTGTTGGAAATTAGTTATGCATATTTGAGAAAATTCTAAAATACTTTGGATAATTGGAGGAGAAGATATTTATAAGGAGCCTTCTAGTGAATATGAGAAGATTGGAATATTCTAGAAGGATTAAGAGGTTGCTTGAAGAGGGTTGCTATTTACCATGGAAAGGTCTTAGAGTAACTTTCAAGAATGGGTGCAAGAATGGGTGTTAGGGTCGTCGCCTTAAAATCTCTACA
>frag0001295_frag0001296
GCATTGAGGGGGAGTGTTGAGAAATAGCAATGCATACTTGAGAATATTCTGGAACACTCTGGATAATTGGAGAGTGTTGGGAATTAGTTATGCATTTTTGAGAAAATTCTAGAATACTTTACTTAATTGGAGGAGAATATATTTACAAGGAGCCTTCTAGTGAATATGAGAAGATTGGAATATTCTAAAAGGCTTGAGAGGTTGTTTGAAGAGGGTTGCTATTTACCATGGAAAGGTCTTAGGGTAACTTTCAAGAGTGGGTGCAAGAGTGGGTGTTAGGGTCGTCGCCCGTAAAATCTCTACA

I'm using the exemplar name as family name, so all hits to a certain exemplar
automatically belong to one family.
"""

import timeForTE
import os
from itertools import groupby

def isLTR_RT(fraglist):
	if all(frag.LTRfeature == 'LTR' or frag.LTRfeature == 'int' for frag in fraglist):
		return(True)
	else:
		return(False)

def has_multiple_LTRs(fraglist):
	if sum(frag.LTRfeature == 'LTR' for frag in fraglist) > 1:
		return(True)
	else:
		return(False)

def int_present(fraglist):
	if any(frag.LTRfeature == 'int' for frag in fraglist):
		return(True)
	else:
		return(False)

def LTR_present(fraglist):
	if any(frag.LTRfeature == 'LTR' for frag in fraglist):
		return(True)
	else:
		return(False)


def concatenate_if_need_be(myfraglist): #returns a list like ['frag0001292_frag0001293', 'frag0001295_frag0001296']
	LTRorINT = [frag.LTRfeature == 'LTR' for frag in myfraglist]
	grouped_L = [(k, sum(1 for i in g)) for k,g in groupby(LTRorINT)]
	n5p = []
	n3p = []
	to_append = []
	if len(grouped_L) == 3:
		n5p = grouped_L[0][1]
		n3p = grouped_L[2][1]
	if len(grouped_L) == 2 and grouped_L[0][0] == False:
		n5p = 0
		n3p = grouped_L[1][1]
	if len(grouped_L) == 2 and grouped_L[0][0] == True:
		n5p = grouped_L[0][1]
		n3p = 0
	if n5p == 1:
		to_append.append(myfraglist[0].fragcode)
	if n5p > 1:
		to_append.append( '_'.join([f.fragcode for f in myfraglist[ : n5p ]]) )
	if n3p == 1 and n5p > 0:
		to_append.append( myfraglist[ (grouped_L[0][1] + grouped_L[1][1]) ].fragcode )
	if n3p > 1 and n5p > 0:
		to_append.append( '_'.join([f.fragcode for f in myfraglist[ (grouped_L[0][1] + grouped_L[1][1]) : len(LTRorINT) ]]) )
	if n3p == 1 and n5p == 0:
		to_append.append( myfraglist[ grouped_L[0][1] ].fragcode )
	if n3p > 1 and n5p == 0:
		to_append.append( '_'.join( [f.fragcode for f in myfraglist[ grouped_L[0][1] : ] ] )  )
	return(to_append)


hella_data_structures = timeForTE.read_fragments_classified()
TEcopyobjs = hella_data_structures[0] #indexed by TEcopy code
fraglists = hella_data_structures[6] #indexed by TEcopy code, a list of TEfrag objects

TEsByFamily = {}  #keys are RMnames (e.g. 'Gypsy-57_BD#LTR/Gypsy', 'Copia-36_BD#LTR/Copia') and value is a list of TEcodes 
for TEcode, TEobj in TEcopyobjs.items():
	if TEobj.RMname not in TEsByFamily:
		TEsByFamily[TEobj.RMname] = [TEcode]
	elif TEobj.RMname in TEsByFamily:
		TEsByFamily[TEobj.RMname].append(TEcode)

fragseqs = {} #looks like: {'frag0000001':'TAAACCCTAAACCCT...', 'frag0000002':'GATATTTACAAAG...', ...}
fragseqsfile = open('allfragments.fasta', 'r') 
read = fragseqsfile.readlines()
fragseqsfile.close()
for n in range(len(read)): 
	if read[n].startswith('>'):
		fragseqs[read[n].strip('>').split('::')[0]] = read[n+1].strip('\n')


if not os.path.exists('LTRs_by_family'):  #Won't raise an error, WILL overwrite your files
	os.mkdir('LTRs_by_family')


slurmkey = open('LTRs_by_family/slurm_aliases.txt', 'w')
slurmcount = 1
for RMn, TEcodelist in TEsByFamily.items():
	familydata = {}
	for TEcode in TEcodelist:
		to_keep = []
		myfraglist = fraglists[TEcode]
		if isLTR_RT(myfraglist):
			if LTR_present(myfraglist):
				if int_present(myfraglist):
					result = concatenate_if_need_be(myfraglist)
					for e in result:
						to_keep.append(e)
				elif not int_present(myfraglist):
					if has_multiple_LTRs(myfraglist):
						to_keep.append( '_'.join([f.fragcode for f in myfraglist]) )
					elif not has_multiple_LTRs(myfraglist):
						to_keep.append(myfraglist[0].fragcode)
		if to_keep:
			for LTRfrag in to_keep:
				familydata[LTRfrag] = ''.join( [fragseqs[myfrag] for myfrag in LTRfrag.split('_')] )
	if familydata:
		familyfilename = 'LTRs_by_family/{}.fa'.format(str(slurmcount))
		slurmkey.write('{}\t{}\n'.format(str(slurmcount), RMn.split('#')[0]))
		slurmcount += 1
		outF = open(familyfilename, 'w')
		for fragcode, seq in familydata.items():
			outF.write('>{}\n'.format(fragcode))
			outF.write('{}\n'.format(seq))
		outF.close()

slurmkey.close()



#write_SBATCH(nodes, ntasks, cpus_per_task, time, mem, qos, mailtype, filename)
timeForTE.write_SBATCH('1', '1', '1', '00:45:00', '3G', timeForTE.qos_partial_node, 'fail', 'LTRs_by_family/runmafft.sh')
with open('LTRs_by_family/runmafft.sh', 'a') as batcharray:
	batcharray.write('#SBATCH --array=1-{}\n'.format(str(slurmcount-1)))
	batcharray.write('#SBATCH --error=align1.%A.%a.err\n')
	batcharray.write('#SBATCH --output=align1%A.%a.out\n')
	batcharray.write('\n')
	batcharray.write('einsi --thread -1 --adjustdirectionaccurately ${SLURM_ARRAY_TASK_ID}.fa > ${SLURM_ARRAY_TASK_ID}.aln\n')
	batcharray.write('sleep 2\n')
	batcharray.write('trimal -in ${SLURM_ARRAY_TASK_ID}.aln -gt 0.8 -out ${SLURM_ARRAY_TASK_ID}.trimmed\n')


if not os.path.exists('LTRs_paired'):
	os.mkdir('LTRs_paired')

slurmkey = open('LTRs_paired/slurm_aliases.txt', 'w')
slurmcount = 1
for TEcode, TEobj in TEcopyobjs.items():
	myfraglist = fraglists[TEcode]
	if timeForTE.is_flLTR_RT(myfraglist):
		outF = open('LTRs_paired/{}.fa'.format(str(slurmcount)), 'w')
		for LTRfrag in concatenate_if_need_be(myfraglist):
			outF.write('>{}\n'.format(LTRfrag))
			outF.write('{}\n'.format(''.join( [fragseqs[myfrag] for myfrag in LTRfrag.split('_')] )))
		outF.close()
		slurmkey.write('{}\t{}\n'.format(str(slurmcount), TEcode))
		slurmcount += 1

slurmkey.close()

#write_SBATCH(nodes, ntasks, cpus_per_task, time, mem, qos, mailtype, filename)
timeForTE.write_SBATCH('1', '1', '1', '00:45:00', '3G', timeForTE.qos_partial_node, 'fail', 'LTRs_paired/runmafft.sh')
with open('LTRs_paired/runmafft.sh', 'a') as batcharray:
	batcharray.write('#SBATCH --array=1-{}\n'.format(str(slurmcount-1)))
	batcharray.write('#SBATCH --error=align1.%A.%a.err\n')
	batcharray.write('#SBATCH --output=align1%A.%a.out\n')
	batcharray.write('\n')
	batcharray.write('einsi --thread -1 --adjustdirectionaccurately ${SLURM_ARRAY_TASK_ID}.fa > ${SLURM_ARRAY_TASK_ID}.aln\n')
	batcharray.write('sleep 2\n')
	batcharray.write('trimal -in ${SLURM_ARRAY_TASK_ID}.aln -gt 0.8 -out ${SLURM_ARRAY_TASK_ID}.trimmed\n')

