# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 11:49:30 2018

@author: virginia
"""
import glob
import sys

mydir = sys.argv[1]

slurm2gene1 = {} 
gene1_to_gene2 = {} #keys are diploid genes, values are hybridum genes (or gene1 and gene2)
with open('{}/slurm_aliases.txt'.format(mydir), 'r') as slurmfile:
    for line in slurmfile.readlines():
        slurmID, gene1, gene2 = line.split('\t')[0], line.split('\t')[1], line.split('\t')[2].strip('\n')
        slurm2gene1[slurmID] = gene1
        gene1_to_gene2[gene1] = gene2

fileList = glob.glob('{}/candidatePseudogeneRegions/*.yn'.format(mydir))
outF = open('{}/dNdSbyGene.txt'.format(mydir), 'w')
for fname in fileList:
    f = open(fname, 'r')
    read = f.readlines()
    f.close()
    if len(read) > 0:
        slurmID = fname.split('candidatePseudogeneRegions/')[1].strip('.yn')
        info = list(filter(lambda a: a != '', read[90].split(' '))) #parse the 90th line, which has the stuff we want
        omega = info[6]
        outF.write('{}\t{}\t{}\n'.format(slurm2gene1[slurmID], gene1_to_gene2[slurm2gene1[slurmID]], omega))
outF.close()

