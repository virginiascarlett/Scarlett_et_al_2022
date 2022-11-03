"""
Clean up output from MACSE for PAML/yn00: Censor stop codons and replace '!' with '?'.
"""
import os
import sys

class Alignment(object):
    def __init__(self, gene1, seq1, gene2, seq2):
        self.gene1 = gene1
        self.seq1 = seq1
        self.gene2 = gene2
        self.seq2 = seq2


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




fname = sys.argv[1] #Should look like e.g. candidatePseudogeneRegions/1.cds_macse_NT.fasta 
alias = int(''.join(filter(str.isdigit, fname.split('/')[-1])))

res = read_fasta(fname)
data, geneOrder = res[0], res[1]
myaln = Alignment(geneOrder[0], data[geneOrder[0]], geneOrder[1], data[geneOrder[1]])
res1 = censor_stop_codons(myaln.seq1)
res2 = censor_stop_codons(myaln.seq2)

censored_seq1, PTCs1 = res1[0], res1[1]
censored_seq2, PTCs2 = res2[0], res2[1]

outFprefix = "/".join(fname.split('/')[:-1])

with open('{}/{}.readyForPAML.fasta'.format(outFprefix, alias), 'w') as outF:
    outF.write('>{}\n'.format(myaln.gene1))
    outF.write('{}\n'.format(censored_seq1.replace("!", "?"))) #MACSE gives frameshifts as exclamation points but yn00 can't handle that
    outF.write('>{}\n'.format(myaln.gene2))
    outF.write('{}\n'.format(censored_seq2.replace("!", "?")))

