"""
A script to concatenate RepeatMasker GFFs from individual chromosomes.
In addition to concatenating files, we are changing the third column from
"similarity" (not very informative) to the superfamily of that TE.

This is basically a modified copy of classify_fragments.py that classifies 
everything in the library, not just everything in allfragments.txt.

Does not require any command line arguments, but you must put timeForTE.py 
in your PYTHONPATH like so: 
export PYTHONPATH="${PYTHONPATH}:/global/projectb/scratch/vstartag/TEannotation/annot_Bd21_w_master"
"""


import timeForTE
import glob
import re

def atoi(text):
	return int(text) if text.isdigit() else text

def natural_keys(text):
	return [ atoi(c) for c in re.split(r'(\d+)', text) ]

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


libdict = {}
libname = glob.glob('{}/step4/*lib*.fa'.format(timeForTE.topdir))[0]
TElib = open(libname, 'r')
with open(libname, 'r') as TElib:
	for line in TElib.readlines():
		if line.startswith('>'):
			RMname = line.split()[0].strip('>')
			libdict[RMname] = timeForTE.Exemplar(RMname)

TEnamesdict = {k.split('#')[0]:v for (k, v) in libdict.items()}

alloutfiles = glob.glob('{}/step4/*.gff'.format(timeForTE.topdir))
alloutfiles.sort(key=natural_keys)
finallines = []
extras = []
for f in alloutfiles:
	with open(f, 'r') as gff:
		read = gff.readlines()
		for n in range(3, len(read)):
			linedata = read[n].split('\t')
			if "-rich" in linedata[8]:
				pass
			#else:
			if ":(" and ")n" in linedata[8]:
					#linedata[2] = "low_complexity"
					#finallines.append(linedata)
				pass
			else:
				TEname = linedata[8].replace('Target "Motif:', '').split()[0].strip('\"')
				if TEname in TEnamesdict:
					linedata[2] = classify_exemplar(TEnamesdict[TEname])
					linedata[8] = "ID={}\n".format(TEname)
					finallines.append(linedata)
				# else: #tried to rescue a couple of families that were not in TEnamesdict due to funny name parsing issues... gave up, not worth the effort
				# 	if not TEname in extras: 
				# 		classif = 'could_not_classify'
				# 		# for k in libdict.keys():
				# 		# 	if TEname.rstrip('int').rstrip('-') in k:
				# 		# 		classif=libdict[k].RMname.split('#')[1] 
				# 		linedata[2] = classif
				# 		finallines.append(linedata)
				# 		# extras.append(TEname)


with open('{}_RepeatMasker_{}.gff'.format(timeForTE.genome, libname.split('.')[0].split('/')[-1]), 'w') as outF:
	for line in finallines:
		outF.write("\t".join(line))

