"""
Get some basic summary statistics for our genome's TE annotation.
For each classification category, how many copies are there?
How many full-length LTR-RTs?
Which chromosome has the highest/lowest TE density?

Updated Dec. 2020: calculating genome space covered by a given 
TE type by merging overlaps, instead of just summing the TE lengths.

Updated April 27 2022: Total LTR-RTs and total DNA TEs now reflects
total copies, not total fragments.

TO DO: make sure your solo counts are also reflecting copies, not fragments.
"""

import timeForTE
import sys
import os
import operator

def is_solo(fraglist):
	if all(frag.LTRfeature == 'LTR' for frag in fraglist):
		return(True)
	else:
		return(False)

def sum_frag_lengths(whichdict): #returns a list of values where each value is the length of a TE copy (sum of fragment lengths)
	finalList = []
	for TEcode in whichdict.keys():
		myfraglist = fraglists[TEcode]
		copylen = 0
		for fragobj in myfraglist:
			copylen += int(fragobj.end)-int(fragobj.start)
		finalList.append(copylen)
	return(finalList)


def write_output(outFname, copynumdict, totalbpdict, lengthdict, mychroms, solodict, fldict):
	outF = open(outFname, 'w')
	outF.write('Classification\tNumber_of_Copies\tTotal_Genome_Space(Mb)\tAvg_length\tChr_w_highest_copynum_per_Mb\tChr_w_lowest_copynum_per_Mb\tChrs_avg_copynum_per_Mb\n')
	DNAsum = 0
	LTRsum = 0
	for classif in classificationsOrdered:
		total_copies = sum(copynumdict[classif].values())
		totalMb = round(sum(totalbpdict[classif].values())/1000000, 2)
		avg_len = round( sum(lengthdict[classif].values())/len(lengthdict[classif]) ) #yes it's an average of averages which is not great but whatever.
		copynum_per_Mb = {} #Looks like e.g. {'Bd1': 54.49468237266198, 'Bd2': 57.00942363574175, 'Bd3': 57.64573510007395, 'Bd4': 57.14592154476147, 'Bd5': 57.52679624015757}
		for chrom in mychroms:
			copynum_per_Mb[chrom] = copynumdict[classif][chrom]/chrom_sizes[chrom]*1000000
		chr_w_highest_copynum = '{}:{}'.format(
			max(copynum_per_Mb.items(), key=operator.itemgetter(1))[0], 
			round(max(copynum_per_Mb.items(), key=operator.itemgetter(1))[1], 2))
		chr_w_lowest_copynum = '{}:{}'.format(
			min(copynum_per_Mb.items(), key=operator.itemgetter(1))[0],
			round(min(copynum_per_Mb.items(), key=operator.itemgetter(1))[1], 2))
		chrs_avg_copynum = round( sum(copynum_per_Mb.values())/len(copynum_per_Mb), 2 )
		outF.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(classif, total_copies, totalMb, avg_len, chr_w_highest_copynum, chr_w_lowest_copynum, chrs_avg_copynum) )
		if classif.split('_')[0] == 'DNA':
			DNAsum += int(total_copies)
		if classif.split('_')[0] == 'LTR':
			LTRsum += int(total_copies)
	outF.write('Total_DNA_transposons: {}\n'.format(str(DNAsum)))
	outF.write('Total_LTR_retrotransposons: {}\n'.format(str(LTRsum))) #I THINK THIS IS COUNTING FRAGMENTS AS COPIES, NEED TO FIX
	outF.write('Total_solo_LTRs: {}\n'.format(str(len(solodict))))
	outF.write('Total_full-length_LTR_Retrotransposons: {}\n'.format(str(len(fldict))))
	sololens = sum_frag_lengths(solodict)
	flLTRRTlens = sum_frag_lengths(fldict)
	outF.write('Avg_length_of_solo_LTRs: {}\n'.format( str(round(sum(sololens)/len(sololens), 2)) ))
	outF.write('Avg_length_of_full-length_LTR_Retrotransposons: {}\n'.format( str(round(sum(flLTRRTlens)/len(flLTRRTlens), 2)) ))
	DNA_list_of_chr_dicts = [total for (classif, total) in totalbpdict.items() if 'DNA_transposon' in classif]
	RT_list_of_chr_dicts = [total for (classif, total) in totalbpdict.items() if 'retrotransposon' in classif]
	DNA_list_of_lists = [list(d.values()) for d in DNA_list_of_chr_dicts]
	RT_list_of_lists = [list(d.values()) for d in RT_list_of_chr_dicts]
	DNA_totalbp = sum([sum(L) for L in DNA_list_of_lists])
	RT_totalbp = sum([sum(L) for L in RT_list_of_lists])
	genomesize = sum([v for (k,v) in chrom_sizes.items() if k in mychroms])
	outF.write('Percent_of_genome_space_occupied_by_DNA_transposons: {}\n'.format( str( round( DNA_totalbp / genomesize *100, 4 ) ) ))
	outF.write('Percent_of_genome_space_occupied_by_retrotransposons: {}\n'.format( str( round( RT_totalbp / genomesize *100, 4) ) ))
	outF.close()


classificationsOrdered = [
'LTR_retrotransposon_RLC', 
'LTR_retrotransposon_RLG', 
'LTR_retrotransposon_RLX', 
'non_LTR_retrotransposon',
'DNA_transposon_CACTA',
'DNA_transposon_hAT',
'DNA_transposon_Mutator',
'DNA_transposon_Mariner',
'DNA_transposon_PIF-Harbinger',
'DNA_transposon_MULE-MuDR',
'DNA_transposon_CMC-EnSpm',
'DNA_transposon_MITE-Stowaway',
'DNA_transposon_MITE-Tourist',
'DNA_transposon_MITE-Unclassified',
'DNA_transposon_RC',
'DNA_transposon_Unclassified',
'could_not_classify'
]

hella_data_structures = timeForTE.read_fragments_classified()
TEcopyobjs = hella_data_structures[0]  #indexed by TEcopy code
fragobjs = hella_data_structures[1] #key is fragcode and value is fragobj
exobjs = hella_data_structures[2] #key is TEcopy code, value is a list: [exemplarobj, classification]
fragstocopies = hella_data_structures[3]  #key is fragcode and value is TEcopycode
classifications = hella_data_structures[4] #key is classification, value is a list of frag objs
flLTRRTs = hella_data_structures[5] #indexed by TEcopy code, value is TEcopy object
fraglists = hella_data_structures[6] #indexed by TEcopy code, a list of TEfrag objects
TEcopiesOrdered = hella_data_structures[7]
chroms = hella_data_structures[8]
starts = hella_data_structures[9]
ends = hella_data_structures[10]
names = hella_data_structures[11]




soloLTRs = {} #indexed by TEcopycode, value is a TEcopy object
for onecode, fraglist in fraglists.items():
	if is_solo(fraglist):
		soloLTRs[onecode] = TEcopyobjs[onecode]

soloLTRs_D = {k:v for (k,v) in soloLTRs.items() if v.chrom.startswith('BhD')}
soloLTRs_S = {k:v for (k,v) in soloLTRs.items() if v.chrom.startswith('BhS')}
flLTRRTs_D = {k:v for (k,v) in flLTRRTs.items() if v.chrom.startswith('BhD')}
flLTRRTs_S = {k:v for (k,v) in flLTRRTs.items() if v.chrom.startswith('BhS')}

chrom_sizes = timeForTE.get_chrom_sizes()

#For each of these dicts, key is a classification, value is a dictionary.
#For the subdictionary: key is a chromosome, value is the <stat> of that TE type on that chromosome
#where <stat> is: 
copynum_by_chrom = {} #copy number
length_by_chrom = {} #average length
totalbp_by_chrom = {} #TE bp, overlaps are counted only once
#Looks like e.g.:
#>>> copynum_by_chrom['LTR_retrotransposon_RLC']
#{'Bd1': 4091, 'Bd2': 3371, 'Bd3': 3438, 'Bd4': 2777, 'Bd5': 1647}
#>>> length_by_chrom['LTR_retrotransposon_RLC']
#{'Bd1': 826.0537765827426, 'Bd2': 795.311480272916, 'Bd3': 839.1239092495637, 'Bd4': 884.1256751890529, 'Bd5': 969.4960534304796}
#>>> totalbp_by_chrom['LTR_retrotransposon_RLC']
#{'BhD1': 3075476, 'BhD2': 2573132, 'BhD3': 2795348, 'BhD4': 2298701, 'BhD5': 1422826, 'BhS1': 965741, 'BhS2': 819196, 'BhS3': 1169458, 'BhS4': 800506, 'BhS5': 905435, 'BhS6': 1489774, 'BhS7': 1095827, 'BhS8': 972371, 'BhS9': 821254, 'BhS10': 1157601}


for classif in classificationsOrdered:
	copynum_by_chrom[classif] = {}
	length_by_chrom[classif] = {}
	totalbp_by_chrom[classif] = {}
	for chromosome in timeForTE.chromosomes:
		intervals = [] #looks like e.g. [['3957', '5279'], ['19484', '19627'], ['97577', '97598'], ['174879', '175009'], ['202935', '203022']...]
		#each entry^ represents one fragment
		lengths = [] #[1322, 143, 21, 130, 87]
		#each entry^ represents one TE copy
		for onecode in TEcopiesOrdered:
			if TEcopyobjs[onecode].chrom == chromosome:
				if exobjs[onecode][1] == classif:
					fraglist = fraglists[onecode]
					length = 0
					for fragobj in fraglist:
						intervals.append([fragobj.start, fragobj.end])
						length += int(fragobj.end)-int(fragobj.start)
					lengths.append(length)
		if intervals:
			merged_intervals = timeForTE.merge_intervals(intervals)
			totalbp_by_chrom[classif][chromosome] = sum([ int(pair[1])-int(pair[0]) for pair in intervals ])
		else:
			totalbp_by_chrom[classif][chromosome] = 0
		if lengths:
			copynum_by_chrom[classif][chromosome] = len(lengths)
			length_by_chrom[classif][chromosome] = sum(lengths)/len(lengths)
		else:
			copynum_by_chrom[classif][chromosome] = 0
			length_by_chrom[classif][chromosome] = 0
		




copynum_by_chrom_D = {}
length_by_chrom_D = {}
totalbp_by_chrom_D = {}
copynum_by_chrom_S = {}
length_by_chrom_S = {}
totalbp_by_chrom_S = {}
for classif in classificationsOrdered:
	copynum_by_chrom_D[classif] = {k:v for (k,v) in copynum_by_chrom[classif].items() if k.startswith('BhD')}
	length_by_chrom_D[classif] = {k:v for (k,v) in length_by_chrom[classif].items() if k.startswith('BhD')}
	totalbp_by_chrom_D[classif] = {k:v for (k,v) in totalbp_by_chrom[classif].items() if k.startswith('BhD')}
	copynum_by_chrom_S[classif] = {k:v for (k,v) in copynum_by_chrom[classif].items() if k.startswith('BhS')}
	length_by_chrom_S[classif] = {k:v for (k,v) in length_by_chrom[classif].items() if k.startswith('BhS')}
	totalbp_by_chrom_S[classif] = {k:v for (k,v) in totalbp_by_chrom[classif].items() if k.startswith('BhS')}





write_output('TEstats.txt', copynum_by_chrom, totalbp_by_chrom, length_by_chrom, timeForTE.chromosomes, soloLTRs, flLTRRTs)

if len(timeForTE.chromosomes) == 15:
	write_output('TEstats.BhD.txt', copynum_by_chrom_D, totalbp_by_chrom_D, length_by_chrom_D, [c for c in timeForTE.chromosomes if c.startswith('BhD')], soloLTRs_D, flLTRRTs_D)
	write_output('TEstats.BhS.txt', copynum_by_chrom_S, totalbp_by_chrom_S, length_by_chrom_S, [c for c in timeForTE.chromosomes if c.startswith('BhS')], soloLTRs_S, flLTRRTs_S)




outF2 = open('bp_bychrom.txt', 'w')
outF2.write('Classification\t{}\n'.format('\t'.join(timeForTE.chromosomes)))
for classif, chromdic in totalbp_by_chrom.items():
	outF2.write('{}\t{}\n'.format(classif, '\t'.join( [str(chromdic[c]) for c in timeForTE.chromosomes] )))

outF2.close()




