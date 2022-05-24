"""

STARTED ON THIS BUT NEVER FINISHED, THIS SCRIPT DOES NOT WORK



This script is adapted from my script classify_fragments.py.

Some articles I used for reference:
https://www.nature.com/articles/nrg2165/figures/1
https://www-sciencedirect-com.libproxy.berkeley.edu/science/article/pii/S1055790315000664
https://tehub.org/classification/wicker

see https://tehub.org/classification/wicker
#there are probably more categories I am missing. I just add new ones to the list as I find them.
abbrevs <- list('LTR/unknown'='RLX',
                'LTR/Cassandra' = 'RLX',
                'LTR'='RLX', 
                'LTR/Caulimovirus'='RLX', 
                'LTR/Gypsy'='RLG', 
                'LTR/Gypsy?'='RLG', 
                'LTR/Copia'='RLC', 
                'LTR/Copia?'='RLC', 
                'LTR/Halcyon'='RLA',
                'LTR/Pao'='RLB', 
                'LTR/Retrovirus'='RLR',
                'LTR/ERV1'='RLE',
                'LINE'='RIX',
                'LINE/LINE'='RIX',
                'LINE?'='RIX',
                'LINE/R2'='RIR',
                'LINE/RTE'='RIT',
                'LINE/RTE-BovB'='RIT',
                'LINE/Jockey'='RIJ',
                'LINE/L1'='RIL',
                'LINE?/L1'='RIL',
                'LINE/I'='RII',
                'SINE'='RSX', 
                'SINE/RTE'='RSX',
                'SINE/R2'='RSX',
                'SINE/tRNA'='RST',
                'SINE/tRNA-RTE'='RST',
                'SINE/7SL'='RSL',
                'SINE/5S'='RSS',
                'LTR/DIRS'='RYD',
                'LTR/PLE'='RPP', 
                'PLE'='RPP',
                'DNA'='DXX', 
                'DNA?'='DXX',
                'DNA/Unknown'='DXX',
                'DNA/CACTA'='DTC', 
                'DNA/CMC-EnSpm'='DTC',
                'DNA/CMC'='DTC',
                'DNA/hAT'='DTA',
                'DNA/hAT-Tip100'='DTA',
                'DNA/hAT-Tip100?'='DTA',
                'DNA/hAT-Ac'='DTA', 
                'DNA/hAT-Tag1'='DTA', 
                'DNA/Mutator'='DTM', 
                'DNA/MULE-MuDR'='DTM', 
                'DNA/MULE-MuDR?'='DTM',
                'DNA/Mariner'='DTT', 
                'DNA/Tc1-Mariner'='DTT',
                'DNA/TcMar-Stowaway'='DTT',
                'DNA/TcMar-Pogo'='DTT',
                'DNA/PIF-Harbinger'='DTH', 
                'DNA/PIF-Harbinger?'='DTH',
                'DNA/Harbinger'='DTH', 
                'DNA/Helitron'='DHH',
                'DNA/RC'='DHH',
                'DNA/Crypton'='DYC'
)









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
	Query whether 'PLE' in exemplar or 'Penelope' in exemplar. I didn't find any.


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
		if myexobj.family.startswith('LTRharvest'):
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


def read_fasta(filename): 
	lines = {}
	with open(filename, 'r') as fastafile:
		currentKey = ''
		for line in fastafile.readlines():
			if line.startswith('>'):
				currentKey = line.strip('\n').strip('>')
			else:
				if currentKey in lines:
					lines[currentKey] += line.strip('\n')
				else:
					lines[currentKey] = line.strip('\n')
	return(lines)




fadict = read_fasta('TREP.RM.combined.final')
for header, seq in fadict.items():
	bn = header.split()[0]
	fam = bn.split('#')[0]
	order = bn.split('#')[1].split('/')[0]
	superfam = bn.split('#')[1].split('/')[1]






























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

"""
Split Repbase headers on whitespace then take the stuff after the '#',
here is everything from the monocot library (NEED TO ADD MITE-Tracker, etc.):

LTR/Gypsy?,LTR_retrotransposon_RLG,RLG
LTR/Gypsy,LTR_retrotransposon_RLG,RLG
LTR/Copia,LTR_retrotransposon_RLC,RLC
LTR/Copia?,LTR_retrotransposon_RLC,RLC
LTR,LTR_retrotransposon_RLX,RLX
LTR/Cassandra,LTR_retrotransposon_RLX,RLX
LTR/Halcyon,LTR_retrotransposon_RLX,RLA
LTR/Echo,LTR_retrotransposon_RLX,RLH
LTR/Caulimovirus,non_LTR_retrotransposon,RXX
LINE,non_LTR_retrotransposon,RIX
LINE/LINE,non_LTR_retrotransposon,RIX
LINE?,non_LTR_retrotransposon,RIX
LINE/RTE-BovB,non_LTR_retrotransposon,RIT
Retroposon/L1-derived,non_LTR_retrotransposon,RIL
Retroposon/L1-dep,non_LTR_retrotransposon,RIL
LINE/L1,non_LTR_retrotransposon,RIL
LINE?/L1,non_LTR_retrotransposon,RIL
LINE/R2,non_LTR_retrotransposon,RIR
LINE/RTE,non_LTR_retrotransposon,RIT
LINE/Jockey,non_LTR_retrotransposon,RIJ
LINE/I,non_LTR_retrotransposon,RII
SINE/7SL,non_LTR_retrotransposon,RSL
SINE/5S,non_LTR_retrotransposon,RSS
SINE/tRNA-RTE,non_LTR_retrotransposon,RST
SINE/tRNA-L1,non_LTR_retrotransposon,RST
SINE/tRNA,non_LTR_retrotransposon,RST
SINE,non_LTR_retrotransposon,RSX
SINE/Chronos,non_LTR_retrotransposon,RSX
SINE/R2,non_LTR_retrotransposon,RSX
SINE/Jokey,non_LTR_retrotransposon,RSX
SINE/I,non_LTR_retrotransposon,RSX
SINE/L1,non_LTR_retrotransposon,RSX
SINE/Pan,non_LTR_retrotransposon,RSX
SINE/RTE,non_LTR_retrotransposon,RSX
DNA?,DNA_transposon_Unclassified,DTX
DNA,DNA_transposon_Unclassified,DTX
DNA/CMC-EnSpm,DNA_transposon_Unclassified,DTX
DNA/CMC,DNA_transposon_Unclassified,DTX
DNA/CACTA,DNA_transposon_CACTA,DTC
DNA/Harbinger,DNA_transposon_PIF-Harbinger,DTH
DNA/PIF-Harbinger?,DNA_transposon_PIF-Harbinger,DTH
DNA/PIF-Harbinger,DNA_transposon_PIF-Harbinger,DTH
DNA/Mariner,DNA_transposon_Mariner,DTT
DNA/TcMar-Stowaway,DNA_transposon_MITE-Stowaway,DTT
DNA/Tc1-Mariner,DNA_transposon_Mariner,DTT
DNA/TcMar-Pogo,DNA_transposon_Mariner,DTT
DNA/P,DNA_transposon_P/DTP
DNA/MULE-MuDR,DNA_transposon_Unclassified,DTM
DNA/Mutator,DNA_transposon_Unclassified,DTM
DNA/MULE-MuDR?,DNA_transposon_Unclassified,DTM
DNA/hAT,DNA_transposon_hAT,DTA
DNA/hAT?,DNA_transposon_hAT,DTA
DNA/hAT-Ac,DNA_transposon_hAT,DTA
DNA/hAT-Tip100,DNA_transposon_hAT,DTA
DNA/hAT-Tip100?,DNA_transposon_hAT,DTA
DNA/hAT-Tag1,DNA_transposon_hAT,DTA
DNA/hAT-Charlie,DNA_transposon_hAT,DTA
RC/Helitron,DNA_transposon_RC,DHH
DNA/Helitron,DNA_transposon_RC,DHH
DNA/RC,DNA_transposon_RC,DHH
DNA/Crypton,DNA_transposon_Unclassified,DYC
Satellite/centromeric,could_not_classify,XXX
ARTEFACT,could_not_classify,XXX
snRNA,could_not_classify,XXX
rRNA,could_not_classify,XXX
Unknown,could_not_classify,XXX
Satellite/subtelomeric,could_not_classify,XXX
tRNA,could_not_classify,XXX
Satellite,could_not_classify,XXX

"""

