


import timeForTE


def specific_to_BhD(family, TElist):
	numD, numS, total = sum([TE.chrom.startswith('BhD') for TE in TElist]), sum([TE.chrom.startswith('BhS') for TE in TElist]), len(TElist)
	if total > 1 and numD/total >= 0.90:
		return(True)
	else:
		return(False)

def specific_to_BhS(family, TElist):
	numD, numS, total = sum([TE.chrom.startswith('BhD') for TE in TElist]), sum([TE.chrom.startswith('BhS') for TE in TElist]), len(TElist)
	if total > 1 and numS/total >= 0.90:
		return(True)
	else:
		return(False)


hella_data_structures = timeForTE.read_fragments_classified()

TEcopyobjs = hella_data_structures[0] #indexed by TEcopy code
#fragobjs = hella_data_structures[1] #key is fragcode and value is fragobj

families = {}  #key is RMname and value is a list of TEobjects

for onecode, TEobj in TEcopyobjs.items():
	if TEobj.RMname in families:
		families[TEobj.RMname].append(TEobj)
	elif TEobj.RMname not in families:
		families[TEobj.RMname] = [TEobj]

subg_spec_fams = {} #key is RMname and value is which subgenome it's specific to
subg_spec_num_copies = {} #key is RMname and value is a list: [numD, numS] 
outF = open('subg_specificity_by_family.txt', 'w')
outF.write('famname\tnum_copies_D\tnum_copies_S\n')
for family, TElist in families.items():
	numD, numS, total = sum([TE.chrom.startswith('BhD') for TE in TElist]), sum([TE.chrom.startswith('BhS') for TE in TElist]), len(TElist)
	outF.write('{}\t{}\t{}\n'.format(family, numD, numS))
	# if specific_to_BhD(family, TElist):
	# 	subg_spec_fams[family] = 'BhD'
	# 	subg_spec_num_copies[family] = [numD, numS]
	# 	outF.write('{}\t{}\t{}\n'.format(family, numD, numS))
	# if specific_to_BhS(family, TElist):
	# 	subg_spec_fams[family] = 'BhS'
	# 	subg_spec_num_copies[family] = [numD, numS]
	# 	outF.write('{}\t{}\t{}\n'.format(family, numD, numS))

outF.close()

# KEEfrags = []
# with open('fragments_in_KEEs.bed', 'r') as KEEfile:
# 	KEEfrags = [ line.split('\t')[3].replace('\n', '') for line in KEEfile.readlines() ]

# KEEfams = {} #key is RMname, value is a list of only those TE copy objs from this family that occur in a KEE
# for frag in KEEfrags:
# 	myfragobj = fragobjs[frag]
# 	myTEobj = TEcopyobjs[myfragobj.TEcode]
# 	if myTEobj.RMname in KEEfams:
# 		if myTEobj.onecode in [x.onecode for x in KEEfams[myTEobj.RMname]]:
# 			pass
# 		else:
# 			KEEfams[myTEobj.RMname].append(myTEobj)
# 	else: 
# 		KEEfams[myTEobj.RMname] = [myTEobj]


# trapped = {} #key is TE code, value is TE obj
# outF = open('trappedTEs.txt', 'w')
# outF.write('TEcopycode\tchrom\tstart\tend\tsense\tFamname\tnum_BhD_copies\tnum_BhS_copies\n')
# for RMname, subg in subg_spec_fams.items():
# 	if RMname in KEEfams:
# 		TEobjlist = KEEfams[RMname]
# 		if any(TE.chrom[:3] != subg for TE in TEobjlist): #if there is any KEE TE in this family that is not on the subgenome this family is specific to
# 			for TE in TEobjlist:
# 				if TE.chrom[:3] != subg:
# 					trapped[TE.onecode] = TE
# 					outF.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(TE.onecode, TE.chrom, TE.start, TE.end, TE.sense, TE.RMname, subg_spec_num_copies[TE.RMname][0], subg_spec_num_copies[TE.RMname][1]))

# outF.close()

# outF = open('subg_specific_all.txt', 'w')
# outF.write('TEcode\tchrom\tstart\tend\tfam_specific_to\tfamsize\tfamname\tnotes\n')
# for family, TElist in families.items():
# 	if specific_to_BhD(family, TElist):
# 		for TEobj in TElist:
# 			if TEobj.chrom.startswith('BhD'):
				#outF.write('{}\t{}\t{}\t{}\tBhD\t{}\t{}\ton_orig_subg\n'.format(TEobj.onecode, TEobj.chrom, TEobj.start, TEobj.end, len(TElist), TEobj.RMname))
			#elif TEobj.chrom.startswith('BhS'):
				#if TEobj.onecode in trapped:
				#	outF.write('{}\t{}\t{}\t{}\tBhD\t{}\t{}\ttrapped_invader\n'.format(TEobj.onecode, TEobj.chrom, TEobj.start, TEobj.end, len(TElist), TEobj.RMname))
				#elif TEobj.onecode not in trapped:
				#	outF.write('{}\t{}\t{}\t{}\tBhD\t{}\t{}\tinvader\n'.format(TEobj.onecode, TEobj.chrom, TEobj.start, TEobj.end, len(TElist), TEobj.RMname))
# 	if specific_to_BhS(family, TElist):
# 		for TEobj in TElist:
# 			if TEobj.chrom.startswith('BhS'):
# 				#outF.write('{}\t{}\t{}\t{}\tBhS\t{}\t{}\ton_orig_subg\n'.format(TEobj.onecode, TEobj.chrom, TEobj.start, TEobj.end, len(TElist), TEobj.RMname))
# 			elif TEobj.chrom.startswith('BhD'):
# 				if TEobj.onecode in trapped:
# 					outF.write('{}\t{}\t{}\t{}\tBhS\t{}\t{}\ttrapped_invader\n'.format(TEobj.onecode, TEobj.chrom, TEobj.start, TEobj.end, len(TElist), TEobj.RMname))
# 				elif TEobj.onecode not in trapped:
# 					outF.write('{}\t{}\t{}\t{}\tBhS\t{}\t{}\tinvader\n'.format(TEobj.onecode, TEobj.chrom, TEobj.start, TEobj.end, len(TElist), TEobj.RMname))


# outF.close()

