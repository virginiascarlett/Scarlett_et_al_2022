"""
We need primary transcript lengths to calculate TPMs, but obviously I only have
transcript lengths for the annotated genes, not the pseudogenes. I don't really
see a good solution to this problem, so I will just get primary transcript length
for the pseudogenes and gene length for the pseudogenes. 

This script is a modified version of countsToTPMbasicNEW.py that is just designed
to work on this one experiment where I am tallying counts using a gff that includes
the putative pseudogenes. 
"""

pseudogenelens = {}
#geneorder = []
with open('/global/cscratch1/sd/vstartag/RNAseqBhyb26/Bhyb26gff_w_pseudogenes.gff3', 'r') as gff:
	read = gff.readlines()
	currentgene = ''
	currentgenelen = 0
	for n in range(3, len(read)):
		line = read[n]
		L = line.split('\t')
		if L[2] == 'gene': # and not L[0].startswith('scaffold'):
			geneID = L[8].split('Name=')[1].strip('\n')
			#geneorder.append(geneID)
			if L[1] == 'pseudogene':
				if not currentgene == '':
					pseudogenelens[currentgene] = currentgenelen
				currentgene = geneID
				currentgenelen = 0
		if L[2] == 'CDS' and L[1] == 'pseudogene':
			currentgenelen += int(L[4])-int(L[3])
	pseudogenelens[currentgene] = currentgenelen

#Next, let's make a dictionary that maps the geneID to predicted transcript length
genelens = {}
with open('BhDprimTrlengths.txt', 'r') as BhDlens:
	for line in BhDlens.readlines():
		genelens[line.split('\t')[0]] = line.split('\t')[1].strip('\n')

with open('BhSprimTrlengths.txt', 'r') as BhSlens:
	for line in BhSlens.readlines():
		genelens[line.split('\t')[0]] = line.split('\t')[1].strip('\n')

geneLengthsDict = {**genelens, **pseudogenelens}

#Hmm. Looks like there are 63 genes that are in the gff but not the primary
#transcripts files. Maybe I excluded genes on scaffolds or something. Well I'll just leave those out. 

TPMdicts = [] #a list of dictionaries
for i in range(1,5):
	countsFileName = '{}_HTScounts.pseudos.txt'.format(str(i))
	#Construct a dictionary that maps geneID to count for this library.
	countsDict = {}  #This will be a dictionary with key=geneID and value=count. 
	with open(countsFileName, 'r') as countsFileObj:
		for line in countsFileObj.readlines():
			countsList = line.split('\t')
			if not countsList[0].startswith('__') and not countsList[0].startswith('LT'):   #To accommodate plastid genes and the stats at the end of an HTScounts file, e.g. __no_feature  26675504
				geneID = countsList[0]
				value = countsList[1].strip('\n')
				countsDict[geneID] = value
	#Now we divide the read counts by the length of each transcript to obtain counts per base. Store these in a dictionary.
	cpbDict = {}
	for geneID in countsDict.keys():
		if geneID in geneLengthsDict: #genes not included in geneLengthsDict are not counted. at the moment, this is just genes on scaffolds (i.e. not in BhD or BhS)
			cpb = float(countsDict[geneID]) / float(geneLengthsDict[geneID])
			cpbDict[geneID] = cpb
	#Next, count up all the CPB values in the library and divide this number by 1,000,000. Divide each CPB
	#by this value to obtain TPMs.
	tpmDict = {}
	sumCPBs = sum([cpbDict[geneID] for geneID in cpbDict])
	for geneID in cpbDict.keys():
		geneFraction = cpbDict[geneID] / sumCPBs
		TPM = geneFraction * 1000000.0
		tpmDict[geneID] = round(TPM,3)
	TPMdicts.append(tpmDict)


with open('all_TPMs_w_pseudogenes.txt', 'w') as outF:
	outF.write('Gene\tRoot\tLeaf\tFloret\tCallus\tTotal\n')
	for geneID in geneLengthsDict.keys():
		outF.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(geneID, str(TPMdicts[0][geneID]), str(TPMdicts[1][geneID]), str(TPMdicts[2][geneID]), str(TPMdicts[3][geneID]), str(round( sum([TPMdicts[n][geneID] for n in range(0,4)]), 2 )) ))

