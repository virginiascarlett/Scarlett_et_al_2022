"""
Note that this is meant for gene level counts. My pipeline uses bbtools, a splice-aware aligner, and HTSeq, 
which calculates total gene counts as the sum of exon counts.

Workflow:
CPB = counts per base, or countvalue/genelength for a given gene
sumCPBs = sum of CPBs for all genes in a given library
geneFraction = CPB/sumCPBs for a given gene
TPM = geneFraction x 10^6
TPM = [GOI_count / GOI_length] / [total [GOI_count / GOI_length] for all genes * 1/1,000,000]
    = CPB / sumCPBs * 1,000,000

run with counts file and 'library name' as command line arguments like so:
python3 countsToTPMbasicNEW.py 1_HTScounts.txt 1
"""

import sys

countsFileName = sys.argv[1]
libraryName = sys.argv[2]

#Construct a dictionary that maps geneID to count for each library.
countsDict = {}  #This will be a dictionary with key=geneID and value=count. 
sortedGeneIDs = []
with open(countsFileName, 'r') as countsFileObj:
    for line in countsFileObj.readlines():
        countsList = line.split('\t')
        if not countsList[0].startswith('__') and not countsList[0].startswith('LT'):   #To accommodate plastid genes and the stats at the end of an HTScounts file, e.g. __no_feature  26675504
            geneID = countsList[0]
            value = countsList[1].strip('\n')
            countsDict[geneID] = value
            sortedGeneIDs.append(geneID)

#Next, let's make a dictionary that maps the geneID to predicted transcript length
geneLengthsDict = {}
with open('BhDprimTrlengths.txt', 'r') as BhDlens:
    for line in BhDlens.readlines():
        geneLengthsDict[line.split('\t')[0]] = line.split('\t')[1].strip('\n')
with open('BhSprimTrlengths.txt', 'r') as BhSlens:
    for line in BhSlens.readlines():
        geneLengthsDict[line.split('\t')[0]] = line.split('\t')[1].strip('\n')

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

#Now let's write our TPMs to a tab-delimited output file.
with open('{}_TPMresults.txt'.format(libraryName), 'w') as outputFile:
    for geneID in sortedGeneIDs:
        if geneID in tpmDict: #again, scaffold genes and chloroplast genes might be in the gff/counts but won't get a tpm
            outputFile.write('{}\t{}\n'.format(geneID,tpmDict[geneID]))

