#A script to convert the gff3 output of TransposonPSI to bed format for bedtools' getfasta.
#Meanwhile, it also creates new RepeatMasker-compatible headers.
#Note that if TransposonPSI found multiple fragments belonging to the same "chain" (TE copy),
#I am taking the genomic coordinates of the start of the first fragment 
#and the end of the last fragment, and I'm calling the intervening 
#genomic sequence the exemplar. It looks like the TransposonPSI
#blast hits do cover most of this intervening sequence.
#This is a Frankenstein script made from pieces of other scripts, so it's not well-written

#Run bedtools in such a way that everything comes out 5' to 3'!

#import glob
import timeForTE
import glob

inFile = glob.glob('*.TPSI.allHits.chains.bestPerLocus.gff3')[0]

"""
If you get an error because a weird transposon family appeared, 
(i.e. it's not obvious what order or superfamily it comes from,)
add it to tpsi_to_RM_families with the name as key and value
(e.g. 'ISC1316':'ISC1316'). Google it and if you can tell which
order it's from, go ahead and add it to tpsi_to_RM_orders with 
value = that order, using one of the orders present already 
(e.g. 'LTR' not 'LTR_Retrotransposon', 'ltr', or anything else).
If you can't tell what order it comes from, just use 'Unknown'
as the value of the dictionary entry, or maybe 'DNA?' if you 
think it's a DNA transposon. Finally, add it to 
tpsi_to_RM_superfamilies with value 'Unknown'.

"""

class TPSI_initial_Exemplar(object): 
	def __init__(self, name):
		self.name = name
		self.source = 'TransposonPSI'
		self.seq = ''
		#self.featurecoords = {} #May be empty. If full, looks like e.g. {LTR1:[1, 970], int:[971, 1274], LTR2:[1275, 2243]}
		tpsi_to_RM_families = {'helitronORF':'Helitron', 'gypsy':'Gypsy', 'MuDR_A_B':'MuDR', 'hAT':'hAT', 'cacta':'CACTA', 'TY1_Copia':'Copia', 'LINE':'LINE', 'ltr_Roo':'ltr_Roo', 'ISC1316':'ISC1316', 'mariner':'Mariner', 'ISa':'ISa'}
		L = self.name.split('\t')
		self.chain = L[8].split(';')[0].replace('ID=', '')
		self.keyword = L[8].split(';')[1].replace(' Target=', '')
		if not self.keyword in tpsi_to_RM_families:
			print('ERROR!! {} is an unrecognized keyword. Add it to the tpsi_to_RM_families, _orders, and _superfamilies dictionaries in reformatTransposonPSI.py!'.format(self.keyword))
		self.Family = tpsi_to_RM_families[self.keyword]
	def Order(self):
		ord = ''
		tpsi_to_RM_orders = {'helitronORF':'RC', 'gypsy':'LTR', 'MuDR_A_B':'DNA', 'hAT':'DNA', 'cacta':'DNA', 'TY1_Copia':'LTR', 'LINE':'LINE', 'ltr_Roo':'DNA?', 'ISC1316':'DNA', 'mariner':'DNA', 'ISa':'DNA'}
		if not self.keyword in tpsi_to_RM_orders:
			print('ERROR!! {} is an unrecognized keyword. Add it to the tpsi_to_RM_orders dic in reformatTransposonPSI.py!'.format(self.keyword))
		ord = tpsi_to_RM_orders[self.keyword]
		return(ord)
	def Superfamily(self):
		sf = ''
		tpsi_to_RM_superfamilies = {'helitronORF':'Helitron', 'gypsy':'Gypsy', 'MuDR_A_B':'MULE-MuDR', 'hAT':'hAT', 'cacta':'CACTA', 'TY1_Copia':'Copia', 'LINE':'LINE', 'ltr_Roo':'Unknown', 'ISC1316':'Unknown', 'ISa':'Unknown', 'mariner':'Mariner'} 
		if not self.keyword in tpsi_to_RM_superfamilies:
			print('ERROR!! {} is an unrecognized keyword. Add it to the dictionaries in in the TPSI_initial_Exemplar class in reformatTransposonPSI.py!'.format(self.keyword))
		sf = tpsi_to_RM_superfamilies[self.keyword]
		return(sf)
	def add_to_seq(self, s):
		self.seq += s
	def RMname(self, barcode):  
		RMn = self.name
		if self.Superfamily() == 'Unknown':
			RMn = '{}.{}#{}'.format(self.Family, barcode, self.Order())
		else:
			RMn = '{}.{}#{}/{}'.format(self.Family, barcode, self.Order(), self.Superfamily())
		return(RMn)



exemplars = {}
prelim_copies = {}
with open(inFile, 'r') as f:
	for line in f.readlines():
		L =line.split('\t')
		chainID = L[8].split(';')[0].replace('ID=', '')
		chrom, start, end, sense = L[0], L[3], L[4], L[6]
		if not chainID in exemplars:
			exemplars[chainID] = TPSI_initial_Exemplar(line.strip('\n'))
			prelim_copies[chainID] = [ timeForTE.TEcopy( chainID, chrom, start, end, sense, exemplars[chainID].RMname(chainID) ) ]
		elif chainID in exemplars:
			prelim_copies[chainID].append(timeForTE.TEcopy( chainID, chrom, start, end, sense, exemplars[chainID].RMname(chainID) ))

final_copies = {}
for chainID, copylist in prelim_copies.items():
	copystart, copyend, precopy1 = min([int(p.start) for p in copylist]), max([int(p.end) for p in copylist]), copylist[0]
	final_copies[chainID] = timeForTE.TEcopy( chainID, precopy1.chrom, copystart, copyend, precopy1.sense, exemplars[chainID].RMname(chainID) )

outF = open('TransposonPSIresults.bed', 'w')
for chainID, copyobj in final_copies.items():
	outF.write('{}\t{}\t{}\t{}\t0\t{}\n'.format(
		copyobj.chrom,
		copyobj.start,
		copyobj.end,
		exemplars[chainID].RMname(chainID),
		copyobj.sense
		))

outF.close()










# exemplars = {}
# fragments = {}
# with open(inFile, 'r') as f:
# 	for line in f.readlines():
# 		L =line.split('\t')
# 		chainID = L[8].split(';')[0].replace('ID=', '')
# 		if not chainID in exemplars:
# 			exemplars[chainID] = TPSI_initial_Exemplar(line.strip('\n'), 'TransposonPSI')
# 			fragments[chainID] = [timeForTE.TEfrag(chainID, L[0], L[3], L[4], L[6])]
# 		elif chainID in exemplars:
# 			fragments[chainID].append(timeForTE.TEfrag(chainID, L[0], L[3], L[4], L[6]))

# copies = {}
# for chainID, fraglist in fragments.items():
# 	copystart, copyend, frag1 = min([int(f.start) for f in fraglist]), max([int(f.end) for f in fraglist]), fraglist[0]
# 	copies[chainID] = timeForTE.TEcopy(frag1.barcode, frag1.chrom, copystart, copyend, frag1.sense)

# outF = open('TransposonPSIresults.bed', 'w')
# for chainID, copyobj in copies.items():
# 	outF.write('{}\t{}\t{}\t{}\t0\t{}\n'.format(
# 		copyobj.chrom,
# 		copyobj.start,
# 		copyobj.end,
# 		exemplars[chainID].RMname(chainID),
# 		copyobj.sense
# 		))

# outF.close()















