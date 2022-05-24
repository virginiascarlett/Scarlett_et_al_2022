"""
Split LTR-RT exemplars from TREP to create input files for BLAST.
We will BLAST the two ends against each other to determine internal LTR-int boundaries.

TO DO: ONLY SPLIT EXEMPLARS WITH AN ORIGINAL ENTRY THAT LOOKS LIKE THIS:
>RLG_Tmon_Laura_609E6-1 Triticum monococcum; Retrotransposon, LTR, Gypsy; complete element, non-autonomous; KEY=216
>RLG_Tmon_Romani_AF459088-1 Triticum monococcum; Retrotransposon, LTR, Gypsy; complete element; KEY=238
>RLC_Tmon_TAR1_AF459088-1 Triticum monococcum; Retrotransposon, LTR, Copia; complete element; KEY=241

i.e. has "retrotransposon", "LTR", and "complete element"!

"""
import os
import shutil

if os.path.isdir('TREP_exemplars_to_blast'):
	shutil.rmtree('TREP_exemplars_to_blast')

keyfile = open('TREP.names', 'r')
read = keyfile.readlines()
keyfile.close()

new_to_old = {}
for n in range(1, len(read)):
	linelist = read[n].split('\t')
	new_to_old[linelist[0]] = linelist[1].strip('\n')

os.mkdir('TREP_exemplars_to_blast')
libfile = open('TREP.RMcompatible.preliminary.fasta', 'r')
read = libfile.readlines()
libfile.close()

for n in range(len(read)):
	line = read[n]
	if line.startswith('>'):
		order = line.split('#')[1].split('/')[0]
		if order == 'LTR':
			oldname = new_to_old[line.strip('>').strip('\n')]
			if "complete element" in oldname.split(';')[2]:
				seq = (read[n+1]).strip('\n')
				first_half, second_half = seq[:len(seq)//2], seq[len(seq)//2:] 
				outF1 = open( 'TREP_exemplars_to_blast/{}_first.fa'.format(line.strip('>').split('#')[0]), 'w')
				outF1.write( '{}_firsthalf\n'.format(read[n].strip('\n')) )
				outF1.write( first_half + '\n' )
				outF1.close()
				outF2 = open( 'TREP_exemplars_to_blast/{}_second.fa'.format(line.strip('>').split('#')[0]), 'w')
				outF2.write( '{}_secondhalf\n'.format(read[n].strip('\n')) )
				outF2.write( second_half + '\n' )
				outF2.close()

