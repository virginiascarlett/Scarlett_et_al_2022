"""
Clean up output from MACSE for PAML/yn00.
(Censor stop codons and replace '!' with '?')
"""
import os
import sys

class Alignment(object):
    def __init__(self, gene1, seq1, gene2, seq2):
        self.gene1 = gene1
        self.seq1 = seq1
        self.gene2 = gene2
        self.seq2 = seq2


# def parse_gblocks(filename): 
#     lines = {}
#     genesOrdered = []
#     with open(filename, 'r') as myfile:
#         currentKey = ''
#         for line in myfile.readlines():
#             if line.startswith('>'):
#                 currentKey = line.strip('\n').strip('>')
#                 genesOrdered.append(line.strip('\n').strip('>'))
#             else:
#                 if currentKey in lines:
#                     lines[currentKey] += line.strip('\n').replace(" ", "")
#                 else:
#                     lines[currentKey] = line.strip('\n').replace(" ", "")
#     return(lines, genesOrdered)

def read_fasta(filename): 
    lines = {}
    genesOrdered = []
    with open(filename, 'r') as fastafile:
        currentKey = ''
        for line in fastafile.readlines():
            if line.startswith('>'):
                currentKey = line.split()[0].strip('>').strip('\n')
                genesOrdered.append(currentKey)
            else:
                if currentKey in lines:
                    lines[currentKey] += line.strip('\n')
                else:
                    lines[currentKey] = line.strip('\n')
    return(lines, genesOrdered)


def censor_stop_codons(myseq):
    codons = [myseq[i:i + 3] for i in range(0, len(myseq), 3)]
    anyPTCs = False
    final_seq = []
    for n in range(len(codons)):
        codon = codons[n]
        if codon=="TGA" or codon=="TAA" or codon=="TAG":
            final_seq.append('???')
            if n != len(codons)-1: #if this not is the final codon, then it is a PTC
                anyPTCs = True
        else:
            final_seq.append(codon)
    return( [''.join(final_seq), anyPTCs] )




fname = sys.argv[1] #looks like e.g. candidatePseudogeneRegions/1.cds_macse_NT.fasta OR sample1/candidatePseudogeneRegions/1.cds_macse_NT.fasta OR candidatePseudogeneRegions/1_macse_NT.fasta
slurmID = int(''.join(filter(str.isdigit, fname.split('/')[-1])))

res = read_fasta(fname)
data, geneOrder = res[0], res[1]
myaln = Alignment(geneOrder[0], data[geneOrder[0]], geneOrder[1], data[geneOrder[1]])
res1 = censor_stop_codons(myaln.seq1)
res2 = censor_stop_codons(myaln.seq2)

censored_seq1, PTCs1 = res1[0], res1[1]
censored_seq2, PTCs2 = res2[0], res2[1]

outFprefix = "/".join(fname.split('/')[:-1])

with open('{}/{}.readyForPAML.fasta'.format(outFprefix, slurmID), 'w') as outF:
    outF.write('>{}\n'.format(myaln.gene1))
    outF.write('{}\n'.format(censored_seq1.replace("!", "?"))) #MACSE gives frameshifts as exclamation points but yn00 can't handle that
    outF.write('>{}\n'.format(myaln.gene2))
    outF.write('{}\n'.format(censored_seq2.replace("!", "?")))

#record which alignments had a PTC, mainly for bugtesting purposes
outFname = ''
if len(fname.split('/')) == 3:  #if this is a run of control genes
    outFname = '{}/PTCs_from_macse.txt'.format(fname.split('/')[0])
else:
    outFname = 'PTCs_from_macse.txt'

if PTCs1 or PTCs2:
    if os.path.exists(outFname): 
        with open(outFname, 'a') as PTCfile:
            PTCfile.write('{}\t{}\t{}\n'.format(slurmID, PTCs1, PTCs2) )
    else:
        with open(outFname, 'w') as PTCfile:
            PTCfile.write('alias\tPTC_in_seq1\tPTC_in_seq2\n')
            PTCfile.write('{}\t{}\t{}\n'.format(slurmID, PTCs1, PTCs2) )
    
