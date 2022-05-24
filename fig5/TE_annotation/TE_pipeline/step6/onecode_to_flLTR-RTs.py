"""
A script to parse onecode's {genome}.all.fasta.out_{chr}.elem_sorted.csv file.
The main purpose of this script is to make a unique barcode 
for each fragment (a hit to the TE library) and each assembled TEcopy, 
since these do not have unique names in the onecode output.

outputs are:
(1) allfragments.txt, which stores the barcodes and everything 
else you need to make timeForTE objects. Columns are: (tab-delimited)
fragment code, chromosome, fragment start, fragment end, LTR or int if applicable, TE copy code, TE copy start, TE copy end, RMname
(2) LTRs_from_fl-LTR-RTs.bed, which lists the genomic locations
of the LTR fragments from flLTR-RTs (as determined by onecode),
and whether the LTR if a 5' or 3' LTR. Looks like:
Bd1 1439532	1440005	onecode000269_frag0000298_five
Bd1	1447299	1447773	onecode000269_frag0000300_three
Bd1	5308156	5308404	onecode001237_frag0001358_five
Bd1	5308392	5308432	onecode001237_frag0001359_five
Bd1	5310280	5310543	onecode001237_frag0001361_three
Bd1	5310531	5310572	onecode001237_frag0001362_three


"""
import timeForTE


def onecode_to_dic(filename):
	linedic = {} #a dictionary where key is the header line and value is a list of all the fragment lines under it.
	#If there is only one fragment for that individual TE, the list is empty.
	TEcopies = []
	with open(filename, 'r') as fobj:
		read = fobj.readlines()
		read = list(filter(lambda a: a != '\n', read)) #remove empty lines
		for n in range(1, len(read)):
			if read[n].startswith('###'):
				TEcopies.append(read[n])
				linedic[read[n]] = []
			else:
				linedic[TEcopies[-1]].append(read[n])
	return(linedic)

def make_TEcopy_obj(onecode_line, code):
	L = onecode_line.split('\t')
	chromosome, start, end, sense, elem, family = L[4], L[5], L[6], L[8], L[9], L[10]
	RMname = '#'.join((elem, family))
	if elem.endswith('-LTR') or elem.endswith('-int') or elem.endswith('-INT'): #three characters after the dash
		RMname = '#'.join((elem[:-4], family))
	if elem.endswith('-I'): #one character after the dash
		RMname = '#'.join((elem[:-2], family))
	return( timeForTE.TEcopy(code, chromosome, start, end, sense, RMname) )

def make_TEfrag_obj(header_line_list, subheader_line_list, fragcount, TEcopyobj): 
	objlist = []
	tempcount = fragcount
	if subheader_line_list:
		for line in subheader_line_list:
			L = line.split('\t')
			start, end = L[5], L[6]
			feat = 'NA'
			if L[9].endswith('-LTR'):
				feat = 'LTR'
			if L[9].endswith('-int') or L[9].endswith('-I') or L[9].endswith('-INT'):
				feat = 'int'
			objlist.append( timeForTE.TEfrag('frag{}'.format(str(tempcount).zfill(7)), TEcopyobj.chrom, start, end, TEcopyobj.onecode, feat) )
			tempcount += 1
	else:
		feat = 'NA'
		Lh = header_line_list.split('\t')
		if Lh[9].endswith('-LTR'):
			feat = 'LTR'
		if Lh[9].endswith('-int') or Lh[9].endswith('-I') or Lh[9].endswith('-INT'):
			feat = 'int'
		objlist.append( timeForTE.TEfrag('frag{}'.format(str(tempcount).zfill(7)), TEcopyobj.chrom, TEcopyobj.start, TEcopyobj.end, TEcopyobj.onecode, feat) )
		tempcount += 1
	return(objlist, tempcount)



"""
Give each TE copy a unique code which we will use to 
access the TE copy or its fragments in two dictionaries
"""
allTEcopies = {} #keys are onecode codes e.g. 'onecode000001', values are TEcopy objs
fragsbyTEcode = {} #keys are onecode codes e.g. 'onecode000001', values are lists of TEfrag objs
TE_counter = 1
frag_counter = 1
for chrname in timeForTE.chromosomes:
    onecodefile = '{}.all.fasta.out_{}.elem_sorted.csv'.format(timeForTE.genome, chrname)
    linedic = onecode_to_dic(onecodefile) #a dictionary where key is the header line and value is a list of all the fragment lines under it.
	#If there is only one fragment for that individual TE, the list is empty.
    for header, subheads in linedic.items():
        TEcode = 'onecode{}'.format(str(TE_counter).zfill(6))
        TEobj = make_TEcopy_obj(header, TEcode)
        allTEcopies[TEcode] = TEobj
        newlist, newcount = make_TEfrag_obj(header, subheads, frag_counter, TEobj)
        fragsbyTEcode[TEcode] = newlist
        frag_counter = newcount
        TE_counter += 1

"""
Note that we have changed the 'RMname' for LTR-RT exemplars: 
if the RMname has -LTR or -int or -I, we have removed that part. 
So for example:

TElib entries:
>Copia-51_BD-I#LTR/Copia RepbaseID: Copia-51_BD-IXX
>Copia-51_BD-LTR#LTR/Copia RepbaseID: Copia-51_BD-LTRXX

RepeatMasker.out file:
  1521   12.5  1.6  0.4  Bd1        5308156  5308404 (69763141) + Copia-51_BD-LTR         LTR/Copia                   1    252   (325)  3529  
   317    7.3  2.4  0.0  Bd1        5308392  5308432 (69763113) + Copia-51_BD-LTR         LTR/Copia                 536    577     (0)  3529 *
  7883    6.4  2.3  2.0  Bd1        5308402  5310294 (69761251) + Copia-51_BD-I           LTR/Copia                   1   2003     (0)  3530  
  1502   12.7  1.6  0.5  Bd1        5310280  5310543 (69761002) + Copia-51_BD-LTR         LTR/Copia                   1    252   (325)  3530 *
   381    2.4  0.0  0.0  Bd1        5310531  5310572 (69760973) + Copia-51_BD-LTR         LTR/Copia                 536    577     (0)  3530 *
  3436    9.0  5.1  4.0  Bd1        6264531  6265099 (68806446) + Copia-51_BD-LTR         LTR/Copia                   1    575     (2)  4242  
   638   13.6 10.0  0.0  Bd1        7385949  7386058 (67685487) C Copia-51_BD-LTR         LTR/Copia                 (0)    577     457  5113  

onecode output file:
###1521/317/7883/1502/381	7.591	2.123	1.625	Bd1	5308156	5310572	2591	+	Copia-51_BD-LTR	LTR/Copia	NA	NA	NA	3529/3529/3530/3530/3530	5	0.821
1521	12.5	1.6	0.4	Bd1	5308156	5308404	252	+	Copia-51_BD-LTR	LTR/Copia	1	252	(325)	3529	1
317	7.3	2.4	0.0	Bd1	5308392	5308432	42	+	Copia-51_BD-LTR	LTR/Copia	536	577	(0)	3529	1
7883	6.4	2.3	2.0	Bd1	5308402	5310294	2003	+	Copia-51_BD-I	LTR/Copia	1	2003	(0)	3530	1
1502	12.7	1.6	0.5	Bd1	5310280	5310543	252	+	Copia-51_BD-LTR	LTR/Copia	1	252	(325)	3530	1
381	2.4	0.0	0.0	Bd1	5310531	5310572	42	+	Copia-51_BD-LTR	LTR/Copia	536	577	(0)	3530	1
###3436	9.0	5.1	4.0	Bd1	6264531	6265099	575	+	Copia-51_BD-LTR	LTR/Copia	1	575	(2)	4242	1	0.997
###638	13.6	10.0	0.0	Bd1	7385949	7386058	121	C	Copia-51_BD-LTR	LTR/Copia	(0)	577	457	5113	1	0.210

allfragments.txt:
frag0001358	Bd1	5308156	5308404	LTR	onecode001237	5308156	5310572	Copia-51_BD#LTR/Copia
frag0001359	Bd1	5308392	5308432	LTR	onecode001237	5308156	5310572	Copia-51_BD#LTR/Copia
frag0001360	Bd1	5308402	5310294	int	onecode001237	5308156	5310572	Copia-51_BD#LTR/Copia
frag0001361	Bd1	5310280	5310543	LTR	onecode001237	5308156	5310572	Copia-51_BD#LTR/Copia
frag0001362	Bd1	5310531	5310572	LTR	onecode001237	5308156	5310572	Copia-51_BD#LTR/Copia

"""
outF = open('allfragments.txt', 'w')
outF.write('fragcode\tchrom\tfrag_start\tfrag_end\tLTRfeature\tTEcopycode\tTEcopy_start\tTEcopy_end\tstrand\texemplar\n')
for code, fragobjlist in fragsbyTEcode.items():
	for fragobj in fragobjlist:
		outF.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
			fragobj.fragcode,
			fragobj.chrom,
			fragobj.start,
			fragobj.end,
			fragobj.LTRfeature,
			code,
			allTEcopies[code].start,
			allTEcopies[code].end,
			allTEcopies[code].sense,
			allTEcopies[code].RMname)
		)

outF.close()


outF = open('LTRs_from_fl-LTR-RTs.bed', 'w')

for code, TEobj in allTEcopies.items():
	if 'LTR' in TEobj.RMname:
		if fragsbyTEcode[code]:
			if timeForTE.is_flLTR_RT(fragsbyTEcode[code]):
				are_we_there_yet = 'no'  #tells us whether we have hit the int fragment(s) yet
				for fragobj in fragsbyTEcode[code]:
					if fragobj.LTRfeature == 'LTR':
						if are_we_there_yet == 'no':
							outF.write('{}\t{}\t{}\t{}_{}_five\n'.format(fragobj.chrom, fragobj.start, fragobj.end, code, fragobj.fragcode))
						elif are_we_there_yet == 'yes':
							outF.write('{}\t{}\t{}\t{}_{}_three\n'.format(fragobj.chrom, fragobj.start, fragobj.end, code, fragobj.fragcode))
					elif fragobj.LTRfeature == 'int':
						are_we_there_yet = 'yes'

outF.close()


"""
Looks like:

Bd1 1439532	1440005	onecode000269_frag0000298_five
Bd1	1447299	1447773	onecode000269_frag0000300_three
Bd1	5308156	5308404	onecode001237_frag0001358_five
Bd1	5308392	5308432	onecode001237_frag0001359_five
Bd1	5310280	5310543	onecode001237_frag0001361_three
Bd1	5310531	5310572	onecode001237_frag0001362_three
Bd1	5908408	5908670	onecode001386_frag0001522_five
Bd1 5914831	5915095	onecode001386_frag0001524_three

"""
