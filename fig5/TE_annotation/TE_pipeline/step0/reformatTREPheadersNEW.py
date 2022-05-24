#OKAY SO NEW PLAN: 
#STEP0 SHOULD LOOK LIKE THIS: 
#(1) SIMPLIFY TREP HEADER NAMES, STORE IN A FILE. MEANWHILE, CREATE A NEW VERSION OF THE TREP DB WITH THE SIMPLE NAMES (THIS SCRIPT)
#(2) CAT TREP WITH NEW NAMES TO REPEATMASKER LIB
#(3) REMOVE REDUNDANCY 
#New names, which I call "RepeatMasker-compatible", or "RMnames" are unique


#See the bottom of this script for some useful albeit messy notes on RepeatMasker format


class TREP_initial_Exemplar(object): 
	def __init__(self, name):
		self.name = name
		self.seq = ''
		#self.featurecoords = {} #May be empty. If full, looks like e.g. {LTR1:[1, 970], int:[971, 1274], LTR2:[1275, 2243]}
		self.Family = self.name.split('_')[2]
		self.classinfo = self.name.split(';')[1]
		self.SFcode = self.name.split('_')[0]
	def Order(self):
		ord = ''
		if self.classinfo == ' ':
			ord = 'Unknown'
		else:
			ord = self.classinfo.split(',')[1]
			if self.classinfo.split(',')[1] ==' non-LTR (SINE)':
				ord = 'SINE'
			elif self.classinfo.split(',')[1] ==' TIR' or self.classinfo.split(',')[1] ==' Helitron':
				ord='DNA'
			elif self.classinfo.split(',')[1] ==' LTR':
				ord='LTR'
			elif self.classinfo.split(',')[1] ==' unknown':
				ord='Unknown'
		return(ord)
	def Superfamily(self):
		sf = ''
		if self.classinfo == ' ':
			sf = 'Unknown'
		else:
			sf = self.classinfo.split(',')[2].strip(' ')
			if sf == 'unknown':
				sf = 'Unknown'
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

#Create a dictionary that associates each TREP TE with an arbitrary
#ID number which will make the final headers unique.


TEobjs = []
with open('trep-db_nr_Rel-16.nulls-removed.fasta', 'r') as origFasta:
	for line in origFasta.readlines():
		if line.startswith('>'):
			TREPheader = line.strip('\n').strip('>')
			TEobjs.append(TREP_initial_Exemplar(TREPheader))
		else:
			TEobjs[-1].add_to_seq(line.strip('\n'))

keyfile = open('TREP.names', 'w')
newTREP = open('TREP.RMcompatible.preliminary.fasta', 'w')
keyfile.write('RM-compatible_header\toriginal_TREP_header\n')
counter = 1
for e in TEobjs:
	exCode = 'trep{}'.format(str(counter).zfill(5))
	keyfile.write('{}\t{}\n'.format(e.RMname(exCode), e.name))
	newTREP.write('>{}\n{}\n'.format(e.RMname(exCode), e.seq))
	counter += 1

keyfile.close()
newTREP.close()





# Examples of TREP headers:
# >XXX_Taes_unnamed_AF472572-1 Triticum aestivum; unknown, unknown, unknown; tandem repeat, partial unit; KEY=13
# >DTM_Tmon_Pilifon_AF459088-1 Triticum monococcum; DNA-transposon, TIR, Mutator; WARNING: likely TF gene derived from transposase; KEY=244
# >RLC_Hvul_Bianca_AF521177-1 Hordeum vulgare; Retrotransposon, LTR, Copia; complete element; KEY=174
# >DTH_Hvul_Islay_296827-1 Hordeum vulgare; DNA-transposon, TIR, Harbinger; Tourist MITE, complete element; KEY=424
# >DTC_Atau_Mandrake_AF446141-1 Aegilops tauschii; DNA-transposon, TIR, CACTA; complete element; KEY=691
# >RLC_Hvul_Lena_AY853252-1 Hordeum vulgare; Retrotransposon, LTR, Copia; int domain CDS; KEY=1034  #actually this is the only header with "int domain CDS"
# >RLC_Taes_Oref_1020F19-1 Triticum aestivum; Retrotransposon, LTR, Copia; 3' end truncated?; KEY=962
# >RLX_Null_Cassandra_consensus-1 Triticeae species; Retrotransposon, LTR, unknown; consensus sequence; KEY=1712
# >RLG_Cele_GypsyA_RND-1 Caenorhabditis elegans; Retrotransposon, LTR, Gypsy; fragment; KEY=1774
# >RIX_Taes_Adair_3Bsr2-36 Triticum aestivum; Retrotransposon, non-LTR (SINE), unknown; fragment; KEY=2044
# >DTT_Bdis_BdisStowawayT_consensus-1 Brachypodium distachyon; DNA-transposon, TIR, Mariner; consensus sequence; KEY=2121
# >RLC_Hvul_TAR4_AY853252-2 Hordeum vulgare; Retrotransposon, LTR, Copia; solo-LTR?; KEY=1056
# >DHH_Osat_Paddington_consensus-1 Oryza sativa; DNA-transposon, Helitron, Helitron; consensus sequence; KEY=1577

# Examples of RepeatMasker headers:
# >IS150#"ARTEFACT", 
# >(CAAGT)n#Simple_repeat RepbaseID: XX
# >DNA9-78_OS#DNA/MULE-MuDR? RepbaseID: DNA9-78_OSXX
# >CACTA-G#DNA/CMC-EnSpm
# >STOWAWAY24_OS#DNA/TcMar-Stowaway RepbaseID: STOWAWAY24_OSXX
# >BALN1_HV#LINE/L1 RepbaseID: BALN1_HVXX
# >THRIA#DNA/hAT? RepbaseID: THRIAXX
# >TREP43#DNA/hAT-Ac RepbaseID: TREP43XX
# >SZ-11_LTR#LTR RepbaseID: SZ-11LTRXX
# >SZ-23_LTR#LTR/Gypsy RepbaseID: SZ-23LTRXX
# >p-SINE1_OS#SINE RepbaseID: p-SINE1_OSXX
# >Helitron1_OS#RC/Helitron RepbaseID: HELITRON1_OSXX
# >SINE6_OS#SINE/tRNA RepbaseID: SINE6_OSXX
# >Afa_HV#Satellite RepbaseID: AFA_HVXX
# >TREP107#Satellite/subtelomeric RepbaseID: TREP107XX
# >HARB-N6_SBi#DNA/PIF-Harbinger RepbaseID: HARB-N6_SBiXX
# >SZ-10old_LTR#LTR/Copia? RepbaseID: SZ-10LTRXX
# >Copia1_HV-int#LTR/Copia RepbaseID: Copia1_HV_IXX


# To summarize, it looks like TREP format is:
# A_B_C_D Genus species; class, order, superfam; intact-ness status; KEY=1234
# where:
# A = three-letter superfamily code (Wicker classification)
# B = species code
# C = family
# D = some sort of code
# class might be "DNA-transposon", "Retrotransposon" or "unknown"
# order might be "TIR", "LTR", "non-LTR (SINE)", "Helitron", maybe there are others I'm missing
# Superfamily may be "Mutator", "CACTA", "Gypsy", "Copia?", "unknown", "Helitron", "L1", or others
# However, some entries lack any class/order/superfam info, or they may only have order. 
# Here are all the possible results of line.split(';')[1]:
#  DNA-transposon, TIR, Harbinger
#  DNA-transposon, Helitron, Helitron
#  Retrotransposon, LTR, Echo
#  Retrotransposon, LTR, unknown
#  Retrotransposon, non-LTR (SINE), L1
#  Retrotransposon, non-LTR (SINE), unknown
#  LTR
#  DNA-transposon, TIR, Mutator
#  Retrotransposon, non-LTR (SINE), Jokey
#  Retrotransposon, LTR, Gypsy
#  unknown
 
#  Retrotransposon, non-LTR (SINE), Pan
#  Helitron
#  DNA-transposon, TIR, Mariner
#  Retrotransposon, non-LTR (SINE), I
#  DNA-transposon, TIR, CACTA
#  DNA-transposon, TIR, hAT
#  Retrotransposon, LTR, Copia
#  TIR
#  DNA-transposon, TIR, unknown
#  Retrotransposon, non-LTR (SINE), Chronos
#  non-LTR (SINE)
#  unknown, unknown, unknown
#  Retrotransposon, LTR, Halcyon
#  DNA-transposon, unknown, unknown
#  Retrotransposon, non-LTR (SINE), R2

# Meanwhile, RepeatMasker format is:
# >fam#order
# OR
# >fam#order/superfam
# where order may be any of the following:
# LINE?
# LTR
# ARTEFACT
# rRNA
# DNA
# Simple_repeat
# Unknown
# Retroposon
# Satellite
# DNA?
# snRNA
# tRNA
# RC
# LINE
# SINE
# Low_complexity
# And the superfam may be 'Gypsy', 'Copia', 'Caulimovirus', or many other things, and may have a question mark.
# I think if the family name ends with -int or _int that means it's an internal sequence,
# and if it ends with -LTR or _LTR, that means it's just an LTR sequence.



