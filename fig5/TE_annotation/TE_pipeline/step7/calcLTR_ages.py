"""
A script to calculate the "age" (time since insertion) of individual 
fl-LTR-RTs discovered by my pipeline. I am following the instructions
of Wicker et al. 2018 (TEs in bread wheat):
age = divergence/(2*mutation rate) with a mutation rate of 1.3*10^–8

This goes in the dir called LTRs_paired, and just works on EMBOSS
distance matrices produced from aligning the 5' and 3' LTR of a single TE.
Yes the ages are quite strange, with some of my TEs being older than life itself.
(See Sep. 13, 2020 entry in my notebook.)
But maybe I'll still plot these and see if I can see anything.

Wicker 2018 also said:
The lifespan of an individual LTR-RT subfamily was defined as the 
5th to 95th percentile interval between the oldest and youngest insertions.
"""

import glob
import math

outFile = open('finalAges.txt', 'w')
outFile.write('TEcopy\tInsertionTime\n')
mutation_rate = 1.3*pow(10, -8)
famstowrite_l = []
famstowrite_d = {}
for dm in glob.glob('*.distmat'):
	code = dm.split('.')[0]
	f = open(dm, 'r')
	read = f.readlines()
	f.close()
	pID = read[8].split('\t')[2].strip()
	if not 'nan' in pID and pID != '-0.00':
		ID = (float(pID))/100
		age = round(ID/(2*mutation_rate))
		famstowrite_d[code] = age
		famstowrite_l.append(code)
		#outFile.write('{}\t{}\n'.format(code, str(age)))

famstowrite_l = sorted(famstowrite_l, key=lambda x: famstowrite_d[x])
slurm_key = {}
with open('slurm_aliases.txt', 'r') as slurmfile:
	for line in slurmfile.readlines():
		slurm_key[line.split('\t')[0]] = line.split('\t')[1].strip('\n')
for f in famstowrite_l:
	outFile.write('{}\t{}\n'.format(slurm_key[f], str(famstowrite_d[f])))
