
import os
from itertools import groupby, chain

#######################   EDIT HERE   ###########################
topdir = '/global/projectb/scratch/vstartag/TEannotation/Bd21' #the directory that contains step1, step2, etc.
treplibdir = 'global/projectb/scratch/vstartag/TEannotation/Bd21/step0' #The location of the file TREP.RM.combined.final,
#can be anywhere
chromosomes = ['Bd1', 'Bd2', 'Bd3', 'Bd4', 'Bd5'] #I'm not sure what will happen if you include scaffolds
genome = 'Bd21'  #Do not include ANY special characters!! Letters and numbers only!! RepeatModeler will get mad
unmasked_genome_file = 'Bdistachyon_314_v3.0.unmasked.fa' #store this in your top-level driectory,
#the same location where this script is stored
softwaredir = '/global/projectb/sandbox/plant/hybridum/software'
qos_whole_node = 'genepool' #Will run on this qos if it needs 64 cpus
qos_partial_node = 'genepool_shared' #will run on this qos if it needs fewer than 64 cpus
#^These don't have to be different
acct = 'plntanot' #slurm -A
nodetype = 'haswell' #slurm -C
email = 'vstartaglio@lbl.gov'
###############################################################

class TEcopy(object):  
	def __init__(self, onecode, chrom, start, end, sense, RMname):
		self.onecode = onecode
		self.chrom = chrom 
		self.start = start #start in genome
		self.end = end #end in genome
		self.sense = sense if sense=='+' or sense=='-' else '-' #KEEP for reformatTransposonPSI.py! change 'C' to '-'
		self.RMname = RMname


class TEfrag(object):
	def __init__(self, fragcode, chrom, start, end, TEcode, LTRfeature):
		self.fragcode = fragcode  #a code unique to this fragment
		self.chrom = chrom
		self.start = start  #start in genome
		self.end = end  #end in genome
		self.TEcode = TEcode #parent's onecode code
		self.LTRfeature = LTRfeature  # should be 'NA' for any fragment except LTR-RT fragments


#RepeatMasker format is:
# >fam#order
# OR
# >fam#order/superfam
#Superfam may be 'Gypsy','?', 'Copia', 'Caulimovirus', etc....
#Order may be 'LINE', 'LINE?', 'LTR', 'ARTEFACT', 'DNA', 'Retroposon', etc...
#In some cases TransposonPSI only gives family and order, so 
#in those cases I am setting superfam to be whatever order is
class Exemplar(object): 
	def __init__(self, RMname):
		self.RMname = RMname
		self.seq = ''
		self.family = self.RMname.split('#')[0].strip('>')
		self.order = self.RMname.split('#')[1].split('/')[0]
		self.superfamily = self.RMname.split('/')[1].strip('\n') if '/' in self.RMname else ''
	def add_to_seq(self, s):
		self.seq += s





def fasta_to_exemplar_dict(filename):
	fileobj = open(filename, 'r')
	read = fileobj.readlines()
	fileobj.close()
	Exemplars = {} #key is e.g. Abbie.trep00555#LTR/Gypsy, value is an Exemplar object
	ex_names = []
	for n in range(len(read)):
		line = read[n]
		if line.startswith('>'):  # looks like '>Abiba.trep00567#LTR/Gypsy\n'
			name = line.strip('>').strip('\n')   # looks like 'Abiba.trep00567#LTR/Gypsy'
			ex_names.append(name)
			Exemplars[name] = Exemplar(name)
		else:
			Exemplars[ex_names[-1]].add_to_seq(read[n].strip('\n'))
	return(Exemplars)



def reformat_fasta(inFile, outFile, every): #'every' is an integer
	with open(inFile, 'r') as fobj:
		elems = []
		data = []
		current_elem_data = []
		for line in fobj.readlines():
			if line.startswith('>'):
				if elems:
					data.append(''.join(current_elem_data))
				current_elem_data = []
				elems.append(line)
			elif not line.startswith('>'):
				current_elem_data.append(line.rstrip())
		data.append(''.join(current_elem_data))
		outF = open(outFile, 'w')
		for i in range(len(elems)):
			outF.write(elems[i])
			if every == 0:  #Write the whole sequence on one line
				outF.write('{}\n'.format(data[i]))
			else: #Insert line breaks every 'every' characters, e.g. every 50 characters
				seq = data[i]
				outF.write('{}\n'.format('\n'.join(seq[i:i+every] for i in range(0, len(seq), every))))
		outF.close()

def read_fasta(filename): 
	lines = {}
	with open(filename, 'r') as fastafile:
		currentKey = ''
		for line in fastafile.readlines():
			if line.startswith('>'):
				currentKey = line.strip('\n')#.strip('>')
			else:
				if currentKey in lines:
					lines[currentKey] += line.strip('\n')
				else:
					lines[currentKey] = line.strip('\n')
	return(lines)

def sort_fasta(filename): 
	lines = {}
	with open(filename, 'r') as fastafile:
		currentKey = ''
		for line in fastafile.readlines():
			if line.startswith('>'):
				currentKey = line.strip('\n')
			else:
				if currentKey in lines:
					lines[currentKey] += line.strip('\n')
				else:
					lines[currentKey] = line.strip('\n')
	seqlens = {}
	for k, v in lines.items():
		if len(v) in seqlens:
			seqlens[len(v)].append(k)
		else:
			seqlens[len(v)] = [k]
	ordered_lengths = sorted(list(seqlens.keys()), reverse=True)
	outfilename = '.'.join(filename.split('.')[:-1]+['sorted'])+'.'+filename.split('.')[-1]
	outF = open(outfilename, 'w')
	for n in ordered_lengths:
		for header in seqlens[n]:
			outF.write('{}\n'.format(header))
			outF.write('{}\n'.format(lines[header]))
	outF.close()


def genome_to_chromosomes(genomefile, genomename):
	chroms = {}
	chromsOrdered = []
	with open(genomefile, 'r') as genomefileobj:
		currentChrom = ''
		currentSeq = ''
		for line in genomefileobj.readlines():
			if line.startswith('>'):
				chrN = line.strip('>').strip('\n')
				if currentChrom != '':
					chroms[currentChrom] = currentSeq
				chromsOrdered.append(chrN)
				currentChrom = chrN
				currentSeq = ''
			else:
				currentSeq += line
		chroms[currentChrom] = currentSeq
	#for chrom in chromsOrdered:
	for chrom in chromosomes:
		outF = open('{}.{}.fasta'.format(genomename, chrom), 'w')
		outF.write('>{}\n'.format(chrom))
		outF.write(chroms[chrom].strip('\n'))
		outF.close()

def merge_intervals(intervals): #Expects a list of sublists. Each sublist is a pair of ints. The second int in the pair must be larger than the first!
	intervals_int = [[int(pair[0]), int(pair[1])] for pair in intervals]
	intervals_sorted = sorted(intervals_int) #sorts by first value in pair, e.g. [[2, 6], [3, 9], [8, 13], [9, 11], [15, 18], [15, 20]]
	merged_intervals = [intervals_sorted[0]]
	for pair in intervals_sorted[1:]:
		currentstart, currentend = pair[0], pair[1]
		previousstart, previousend = merged_intervals[-1][0], merged_intervals[-1][1]
		if not currentend > currentstart: #Ignore any pairs where the second int is less than or equal to the first
			pass
		if currentstart == previousstart and currentend <= previousend:
			pass
		if currentstart == previousstart and currentend > previousend:
			merged_intervals[-1][1] = currentend
		if currentstart > previousstart and currentend <= previousend:
			pass
		if currentstart <= previousend and currentend > previousend:
			merged_intervals[-1][1] = currentend
		if currentstart > previousend and currentend > previousend:
			merged_intervals.append(pair)
	return(merged_intervals)

def write_SBATCH(nodes, ntasks, cpus_per_task, time, mem, qos, mailtype, filename): #args should all be strings
	outF = open(filename, 'w')
	outF.write('#!/bin/bash\n')
	outF.write('#SBATCH -N {}\n'.format(nodes))
	outF.write('#SBATCH --ntasks={}\n'.format(ntasks))
	outF.write('#SBATCH --cpus-per-task={}\n'.format(cpus_per_task))
	outF.write('#SBATCH --time={}\n'.format(time))
	outF.write('#SBATCH --mem={}\n'.format(mem))
	outF.write('#SBATCH --qos={}\n'.format(qos))
	outF.write('#SBATCH -A {}\n'.format(acct))
	outF.write('#SBATCH -C {}\n'.format(nodetype))
	outF.write('#SBATCH --mail-user={}\n'.format(email))
	outF.write('#SBATCH --mail-type={}\n'.format(mailtype))
	outF.close()


def parse_blast(keyword, alnfile): 
	"""
	Returns the query or subject coordinates of the first blast HSP 
	(alignment) in a blastn output file. Does not record strand info.
	"""
	allcoords = []
	HSPs = 0
	with open(alnfile, 'r') as alnfobj:
		for line in alnfobj.readlines():
			if line.startswith(" Score"):
				HSPs += 1
			if line.startswith(keyword):
				if HSPs == 1:  #collect coords from the first HSP (alignment) only
					linelist = line.split()
					for e in linelist:
						if e.isdigit():
							allcoords.append( int(e) )
	if allcoords:
		return(min(allcoords), max(allcoords))
	else:
		return(None)

def is_flLTR_RT(fraglist):
	result = False
	if all(frag.LTRfeature == 'LTR' or frag.LTRfeature == 'int' for frag in fraglist):
		if fraglist[0].LTRfeature == 'LTR' and fraglist[-1].LTRfeature == 'LTR':
			if any(frag.LTRfeature == 'int' for frag in fraglist):
				result = True
	return(result)

def get_chrom_sizes(): #returns a dict that looks like {'Bd1': 75071545, 'Bd2': 59130575, ...
	chromseqsdic = {}
	chromsizesdic = {}
	currentchrom = ''
	currentseq = ''
	with open('{}/{}'.format(topdir, unmasked_genome_file), 'r') as myfile:
		for line in myfile.readlines():
			if line.startswith('>'):
				if currentchrom:
					chromseqsdic[currentchrom] = currentseq
				currentchrom = line.strip('>').strip('\n')
				currentseq = ''
			else:
				if line.strip('\n'):  #ignore blank lines
					currentseq += line.strip('\n')
	chromseqsdic[currentchrom] = currentseq
	for c in chromosomes:
		chromsizesdic[c] = len(chromseqsdic[c])
	return(chromsizesdic)

def get_bins(binsize, chromsizes): #returns a dict where key is a chromosome and value is a list containing two lists
	allbins = {}
	for c in chromosomes:
		binstarts = []
		binends = []
		chromsize = chromsizes[c]
		for i in range(1, chromsize, binsize):
			binstarts.append(i)
			binends.append(i+binsize-1)
		binends[-1] = chromsize  #the last bin will not be the same size as the other bins
		allbins[c] = [binstarts, binends]
	return(allbins)


def read_fragments_classified(filename = 'allfragments.classified'):
	inFile = open(filename, 'r')
	read = inFile.readlines()
	inFile.close()
	TEcopyobjs = {}  #indexed by TEcopy code
	fragobjs = {} #key is fragcode and value is fragobj
	exobjs = {} #key is TEcopy code, value is a list: [exemplarobj, classification]
	fragstocopies = {}  #key is fragcode and value is TEcopycode
	classifications = {} #key is classification, value is a list of frag objs
	flLTRRTs = {} #indexed by TEcopy code, value is TEcopy object
	fraglists = {} #indexed by TEcopy code, a list of TEfrag objects
	TEcopiesOrdered_wdups = []
	#data for making a new bed file!:
	chroms = []
	starts = []
	ends = []
	names = []
	strands = []
	for n in range(1, len(read)):
		line = read[n]
		L = line.strip('\n').split('\t')
		fragstocopies[L[0]] = L[5]
		myexobj = Exemplar(L[9]) 
		myTEcopyobj = TEcopy(L[5], L[1], L[6], L[7], L[8], L[9])
		TEcopiesOrdered_wdups.append(L[5])
		myfragobj = TEfrag(L[0], L[1], L[2], L[3], L[5], L[4])
		TEcopyobjs[L[5]] = myTEcopyobj
		fragobjs[L[0]] = myfragobj
		classif = L[11]
		exobjs[L[5]] = [myexobj, classif]
		chroms.append(L[1])
		starts.append(L[2])
		ends.append(L[3])
		names.append(L[0])
		strands.append(L[8])
		if not classif in classifications:
			classifications[classif] = [myfragobj]
		else:
			classifications[classif].append(myfragobj)
		if L[10] == 'True':
			flLTRRTs[L[5]] = myTEcopyobj
		if L[5] in fraglists:
			fraglists[L[5]].append(myfragobj)
		else:
			fraglists[L[5]] = [myfragobj]
	seen = set()
	seen_add = seen.add
	TEcopiesOrdered = [x for x in TEcopiesOrdered_wdups if not (x in seen or seen_add(x))]
	return( [
			TEcopyobjs,  #indexed by TEcopy code
			fragobjs, #key is fragcode and value is fragobj
			exobjs, #key is TEcopy code, value is a list: [exemplarobj, classification]
			fragstocopies,  #key is fragcode and value is TEcopycode
			classifications, #key is classification, value is a list of frag objs
			flLTRRTs, #indexed by TEcopy code, value is TEcopy object
			fraglists, #indexed by TEcopy code, a list of TEfrag objects
			TEcopiesOrdered,
			chroms,
			starts,
			ends,
			names,
			strands
			]
		)



if (__name__ == '__main__'):
	print('all clear!')





