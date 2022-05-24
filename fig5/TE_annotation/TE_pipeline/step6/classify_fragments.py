"""
A script to classify fragments from TREP, Repbase, TransposonPSI, etc.
into one more or less consistent system of TE types.
Run after onecode_to_flLTR-RTs.py.
Produces an output file called allfragments.classified, which is the
input for many of my downstream scripts. 

NOTE TO SELF: Honestly this classification is not that great. I mean, it's fine,
it's taken from Vogel et al. 2010 which Tom Wicker was an author on, and it covers the 
major groups, but I cut some corners, e.g. Pao is just 'RLX'. It would be better to use the entire Wicker
classification scheme as in this figure: https://www.nature.com/articles/nrg2165/figures/1
This paper also has an extremely useful figure (fig. 2):
https://www-sciencedirect-com.libproxy.berkeley.edu/science/article/pii/S1055790315000664
See my script plotBigFamsCount.R which has a better classification scheme than this script.

ALSO NEED TO FIX!!! LTR/Caulimovirus and LTR/DIRS are NOT non-LTR-RTs!

How to classify Repbase 'orders' i.e. # >fam#order/superfam
LTR Retrotransposons
	#LTR/Gypsy
	#LTR/Copia
	#LTR/Copia?
	rnd-6_family-13679#LTR/ERV1
	rnd-6_family-3393#LTR/Pao
	LTR2-ZM#LTR/Cassandra (class I, order LTR, superfamily TRIM, family Cassandra), also 'Cassandra.trep00459#LTR'

	Unclassified LTR-Retrotransposons (These are RLX)
		These have one of the following formats:
			LTRharvest00001#LTR/unknown
			familyname.trep01234#LTR 
			familyname#LTR (i.e. a repbase element)
			rnd-6_family-2397#LTR

Non-LTR Retrotransposons
	#LINE?
	rnd-6_family-1577#LINE/L1
	SINE9_OS#SINE/tRNA
	rnd-6_family-5150#LTR/Caulimovirus
	rnd-6_family-370#LTR/DIRS
	Query whether if 'PLE' in exemplar or 'Penelope' in exemplar. I didn't find any.


DNA Transposons
	#DNA
	#DNA?
	DTX includes:
		CACTA: BdisF.trep00653#DNA/CACTA
		Mutator: BdisEhu.trep00750#DNA/Mutator
		hAT: rnd-6_family-13900#DNA/hAT-Tip100; rnd-4_family-3313#DNA/hAT-Ac; /hAT-Tag1
		PIF/Harbinger: rnd-6_family-8144#DNA/PIF-Harbinger
		EnSpm-N1_TA#DNA/CMC-EnSpm
		MuDR-2_BDi#DNA/MULE-MuDR
		Marimom.trep01233#DNA/Mariner
	DXX includes:
		MITEs:
		Classify BdisStowawayS.trep00646#DNA/Mariner as a Stowaway MITE
		Classify HADES_TA#DNA/TcMar-Stowaway as a Stowaway MITE
		Classify BdisTouristM.trep00714#DNA/Harbinger as a Tourist MITE
		Helitrons and other RC transposons:
		AA.trep00786#DNA/Helitron
		rnd-6_family-1399#DNA/RC

Unclassified DNA transposons:
mite43#DNA/Unknown ??
DNA7-1_BD#DNA

"""

import timeForTE

def classify_LTRretrotransposons(myexobj):
	LTR_retrotransposons = ['Gypsy', 'Gypsy?', 'Copia', 'Copia?', 'ERV1', 'Pao', 'Cassandra']
	non_LTR_retrotransposons = ['DIRS', 'Caulimovirus', 'Penelope', 'PLE']
	if myexobj.superfamily:
		if myexobj.superfamily in LTR_retrotransposons:
			if 'Gypsy' in myexobj.superfamily:
				return('LTR_retrotransposon_RLG')
			if 'Copia' in myexobj.superfamily:
				return('LTR_retrotransposon_RLC')
			else:
				return('LTR_retrotransposon_RLX')
		elif not myexobj.superfamily in LTR_retrotransposons:
			if myexobj.superfamily in non_LTR_retrotransposons:
				return('non_LTR_retrotransposon')
		if 'LTRharvest' in myexobj.family: 
			return('LTR_retrotransposon_RLX')
	if not myexobj.superfamily:
			return('LTR_retrotransposon_RLX')
	return('could_not_classify')


def classify_DNAtransposons(myexobj):
	DNA_transposons_DTX = ['CACTA', 'Mutator', 'hAT-Tip100', 'hAT-Ac', 'hAT-Tag1', 'hAT', 'PIF-Harbinger', 'PIF-Harbinger?', 'CMC-EnSpm', 'MULE-MuDR', 'MULE-MuDR?']
	#DTX_dic = {'CACTA':'DTC', 'Mutator':'DTM', 'hAT-Tip100':'DTA', 'hAT-Ac':'DTA', 'hAT-Tag1':'DTA', 'PIF-Harbinger':'DTH', 'Mariner':'DTC'}
	if myexobj.superfamily:
		if myexobj.superfamily in DNA_transposons_DTX:
			if 'hAT' in myexobj.superfamily:
				return('DNA_transposon_hAT')
			else:
				return('DNA_transposon_{}'.format(myexobj.superfamily.strip('?')))
		if not myexobj.superfamily in DNA_transposons_DTX:
			if myexobj.superfamily == 'TcMar-Mariner' or myexobj.superfamily == 'TcMar-Tc1':
				return('DNA_transposon_Mariner')
			if myexobj.superfamily == 'Mariner':
				if not 'Stowaway' in myexobj.RMname:
					return('DNA_transposon_Mariner')
				if 'Stowaway' in myexobj.RMname:
					return('DNA_transposon_MITE-Stowaway')
			if myexobj.superfamily == 'TcMar-Stowaway':
				return('DNA_transposon_MITE-Stowaway')
			if myexobj.superfamily == 'Harbinger':
				if 'Tourist' in myexobj.RMname:
					return('DNA_transposon_MITE-Tourist')
				else:
					return('DNA_transposon_PIF-Harbinger')
			if myexobj.superfamily == 'RC' or myexobj.superfamily == 'Helitron':
				return('DNA_transposon_RC')
			if myexobj.superfamily == 'Unknown' and 'mite' in myexobj.RMname:
				return('DNA_transposon_MITE-Unclassified')
			if myexobj.superfamily == 'IS3EU':
				return('DNA_transposon_Unclassified')
	if not myexobj.superfamily:
		return('DNA_transposon_Unclassified')
	return('could_not_classify')


def classify_exemplar(myexobj):
	result = 'could_not_classify'
	if myexobj.order == 'LTR':
		result = classify_LTRretrotransposons(myexobj)
	if 'LINE' in myexobj.order or 'SINE' in myexobj.order:
		result = 'non_LTR_retrotransposon'
	if 'DNA' in myexobj.order:
		result = classify_DNAtransposons(myexobj)
	return(result)







fragfile = open('allfragments.txt', 'r')
read = fragfile.readlines()
fragfile.close()

fragmentsOrdered = []
TEcopyobjs = {}  #indexed by TEcopy code
fraglists = {} #indexed by TEcopy code, a list of TEfrag objects
fragobjs = {} #key is fragcode and value is fragobj
exobjs = {} #indexed by TEcopy code
fragstocopies = {}  #key is fragcode and value is TEcopycode
for n in range(1, len(read)):
	line = read[n]
	L = line.strip('\n').split('\t')
	fragmentsOrdered.append(L[0])
	fragstocopies[L[0]] = L[5]
	myexobj = timeForTE.Exemplar(L[9]) 
	myTEcopyobj = timeForTE.TEcopy(L[5], L[1], L[6], L[7], L[8], L[9])
	myfragobj = timeForTE.TEfrag(L[0], L[1], L[2], L[3], L[5], L[4])
	exobjs[L[5]] = myexobj
	TEcopyobjs[L[5]] = myTEcopyobj
	fragobjs[L[0]] = myfragobj
	if L[5] in fraglists:
		fraglists[L[5]].append(myfragobj)
	else:
		fraglists[L[5]] = [myfragobj]

outF = open('allfragments.classified', 'w')
outF.write('fragcode\tchrom\tfrag_start\tfrag_end\tLTRfeature\tTEcopycode\tTEcopy_start\tTEcopy_end\tstrand\texemplar\tisFL\tclassification\n')
for fragcode in fragmentsOrdered:
	onecode = fragstocopies[fragcode]
	isFL = str(timeForTE.is_flLTR_RT(fraglists[onecode]))
	classification = classify_exemplar(exobjs[onecode])
	outF.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
		fragcode,
		fragobjs[fragcode].chrom,
		fragobjs[fragcode].start,
		fragobjs[fragcode].end,
		fragobjs[fragcode].LTRfeature,
		onecode,
		TEcopyobjs[onecode].start,
		TEcopyobjs[onecode].end,
		TEcopyobjs[onecode].sense,
		TEcopyobjs[onecode].RMname,
		isFL,
		classification
		)
	)

outF.close()


	



