"""
A script to reformat the fasta file delivered by LTR_retriever.
The headers delivered by LTR_retriever are close to RepeatMasker format,
but the '_int' and '_LTR' need to be changed to '-int' and '-LTR'
for onecode to recognize entries from the same exemplar.
Furthermore, the first part of the header (the part I call the 'family'
in timeForTE.py) must be the same for int and LTR fragments of one
fl-LTR-RT exemplar.

Note only one LTR is present for each flLTR-RT in the fasta file, same as RepBase

run like so: 
python3 reformatLTR_retriever.py Bhybridum_463_v1.0.unmasked.fa.LTRlib.fa Bhybridum_463_v1.0.unmasked.fa.pass.list


Hmm, looks like for ABR113, 386 flLTR-RTs were identified, but only 
217 of those have two entries in the fasta file. Why are a huge proportion of the
fl-LTR-RTs not in the final fasta? Maybe they are redundant to other
sequences that are present in the fasta file... oh well, guess it's fine
"""

import sys

fastafile, keyfile = sys.argv[1], sys.argv[2]

counter = 1

keyfileobj = open(keyfile, 'r')
read = keyfileobj.readlines()
keyfileobj.close()
flLTRRTdata = {}
for n in range(1, len(read)): 
	line = read[n]
	LTRheaderbase = '>' + '..'.join((line.split()[0].split('..')[0], str(int(line.split()[6][3:].split('..')[0]) - 1))) + '_LTR'
	INTheaderbase = '>' + line.split()[0].split(':')[0] + ':' + line.split()[6][3:] + '_INT'
	newbasename = '>LTRharvest{}'.format(str(counter).zfill(5))
	flLTRRTdata[LTRheaderbase] = newbasename + '-LTR'
	flLTRRTdata[INTheaderbase] = newbasename + '-int'
	counter += 1

fastafileobj = open(fastafile, 'r')
current_elem = ''
seq = ''
finalfasta = {}
for line in fastafileobj.readlines():
	if line.startswith('>'):
		if seq:
			finalfasta[current_elem] = seq
			seq = ''
		if line.split('#')[0] in flLTRRTdata:
			current_elem = '#'.join( (flLTRRTdata[line.split('#')[0]], line.split('#')[1]) )
		if not line.split('#')[0] in flLTRRTdata:
			current_elem = '#'.join( ('>LTRharvest{}-{}'.format(str(counter).zfill(5), line.split('#')[0][-3:]), line.split('#')[1]) )
			counter += 1
	else:
		seq += line

fastafileobj.close()

outfileobj = open('LTRharvest.exemplars.RMready', 'w')
for header, seq in finalfasta.items():
	outfileobj.write(header)
	outfileobj.write(seq)

outfileobj.close()


"""
#second number in the 5_TSD column+1 is the first number in the LTR header
#IN:numbers perfectly match the INT header
BhD1:859050..863418     pass    motif:TGCA      TSD:ATGGC       859045..859049  863419..863423  IN:859535..862936       0.9607  ?       unknown NA      1554106
>BhD1:859050..859534_LTR#LTR/unknown
>BhD1:859535..862936_INT#LTR/unknown


BhS2:14227020..14231983 pass    motif:TGCA      TSD:GTAAC       14227015..14227019      14231984..14231988      IN:14227160..14231843   0.9643  -       unknown LTR     1407402
>BhS2:14227020..14227159_LTR#LTR/unknown
>BhS2:14227160..14231843_INT#LTR/unknown


BhD5:13342787..13356648 pass    motif:TGCA      TSD:ACAAC       13342782..13342786      13356649..13356653      IN:13343264..13356171   0.9139  -       Gypsy   LTR     3519052
>BhD5:13342787..13343263_LTR#LTR/Gypsy
>BhD5:13343264..13356171_INT#LTR/Gypsy

BhS10:2687883..2693204  pass    motif:TGCA      TSD:TCATT       2687878..2687882        2693205..2693209        IN:2688298..2692789     0.9831  -       Gypsy   LTR     656157
>BhS10:2687883..2688297_LTR#LTR/Gypsy
>BhS10:2688298..2692789_INT#LTR/Gypsy

"""