"""
One script in the process of determining internal 
LTR boundaries for TREP exemplars, based on blast 
alignments. 

First we ran reformatTREPheadersNEW.py, which gave
each TREP exemplar a unique barcode and a new header
in RepeatMasker format. It produced two files, called
TREP.names and TREP.RMcompatible.prelilminary.fasta.

Next we ran prepare_TREP_for_blast.py, which
split fl_LTR-RTs in half, and put the
two halves in separate files. (We went to all this
trouble because blastn, as far as I can tell, 
requires query and subject to be in different files.)

After that, we ran trepblast.sh, which aligned the 
two halves to produce alignment files with the
extension '.aln'. 

In this script, we are using the aligned sequences 
to get internal LTR boundaries. We assume that the
outer edges of the exemplar are the correct outer 
boundaries of the LTR.

Output of this script is a new TREP db file
that is the same as the original TREP db file but the 
headers have been changed to RepeatMasker format, AND
each fl-LTR-RT has been broken up into 2 entries.

The two entries for a given fl-LTR-RT with alignable LTRs
have the same name and barcode (e.g. "trep000123"), but for 
one entry, the part preceding the '#' ends with '-LTR', 
and that entry just contains the 5' LTR sequence, and for 
the other entry, the part preceding the '#' ends with 
'-int', and that entry just contains the int sequence.

To get the final coords of the boundaries, we take
the length of the original exemplar, and determine which 
coords in _second.fa (the blast 'subject' file, the second
half of the exemplar,) correspond to which coords in the original file.
It's okay if the final LTR sequences are not the same size.

Note some exemplars didn't have alignable regions.
These are excluded from the output file.

"""

import os
import timeForTE

# def get_internal_boundary(keyword, alnfile):
# 	allcoords = []
# 	with open(alnfile, 'r') as fobj:
# 		for line in fobj.readlines():
# 			if line.startswith(keyword):
# 				linelist = line.split()
# 				for e in linelist:
# 					if e.isdigit():
# 						allcoords.append( int(e) )
# 	if keyword == 'Query ':
# 		if allcoords:
# 			return(max(allcoords))
# 		else:
# 			return(None)
# 	elif keyword == 'Sbjct ':
# 		return(min(allcoords))



exobjs = timeForTE.fasta_to_exemplar_dict('TREP.RMcompatible.preliminary.fasta') #key is e.g. Abbie.trep00555#LTR/Gypsy, value is an Exemplar object



"""
Grab the coordinates of the two internal boundaries.
Grab the 5' LTR and internal sequences and store in separate dicts.
"""
LTRs = {}
ints = {}
for name, obj in exobjs.items():
	basename = name.split('#')[0]   #looks like Abiba.trep00567
	alnfile = '{}/TREP_exemplars_to_blast/{}.aln'.format(os.getcwd(), basename)
	if os.path.isfile(alnfile):
		fiveprime_boundaries = timeForTE.parse_blast('Query ', alnfile)
		if fiveprime_boundaries:
			fiveprime_boundary = timeForTE.parse_blast('Query ', alnfile)[1]
			threeprime_boundary = timeForTE.parse_blast('Sbjct ', alnfile)[0] + len(obj.seq)//2 + 1
			fiveprime_seq = obj.seq[:fiveprime_boundary]
			internal_seq = obj.seq[fiveprime_boundary:(threeprime_boundary-1)]
			LTRs[name] = fiveprime_seq
			ints[name] = internal_seq

outF = open('TREP.RMcompatible.fasta', 'w')
for name, obj in exobjs.items():
	if name in LTRs:
		basename, order_superfam = name.split('#')[0], name.split('#')[1] 
		outF.write('>{}-LTR#{}\n'.format(basename, order_superfam))
		outF.write('{}\n'.format(LTRs[name]))
		outF.write('>{}-int#{}\n'.format(basename, order_superfam))
		outF.write('{}\n'.format(ints[name]))
	else:
		outF.write('>{}\n'.format(name))
		outF.write('{}\n'.format(obj.seq))

outF.close()




#e.g. seq = 'abcdefghi'
#>>> first_half, second_half = seq[:len(seq)//2], seq[len(seq)//2:]
#>>> first_half
#'abcd'
#>>> second_half
#'efghi'
#Let's say the alignment starts at 'g': this is '3' in the blast output, but 7th in seq (using 1-based counting)
#( len(seq)//2 + 1 ) = 4 
#4 + 3 = 7