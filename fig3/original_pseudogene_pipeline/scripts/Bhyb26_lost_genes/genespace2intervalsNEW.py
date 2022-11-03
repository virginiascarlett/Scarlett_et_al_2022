"""
A script to follow up on my 'missing' genes in Bhyb26. To recap: 
I took John Lovell's syntenic orthogroups and identified those where all 
but one (sub)genome was represented. In other words, there was at 
least one gene from rice, both ABR113 subgenomes, and both diploids, 
but only one of the two Bhyb26 subgenomes.

This script takes as input John L.'s syntenic orthogroups and
a file 'lost_genes.txt', produced by my script bin_orthogroups.R.
It delivers a bed file of B. hybridum regions. Each region should
contain the putatively missing gene that is orthologous to the
diploid gene listed in the 'names' column.

The workflow is as follows: 
First, syntenic orthogroups were filtered to discard any that orthogroups
that had more than one gene from the diploid corresponding to the missing 
subgenome. (As it happens, none of them had this issue.) 
Once I had identified the diploid gene corresponding to the 'missing' Bhyb26 
gene, I 'walked' outward along the chromosome in both directions, checking 
whether the each nearby gene had a single ortholog in the appropriate Bhyb26 
subgenome. If it had no orthologs or many orthologs, I skipped it and moved 
on to the next-closest gene. I proceeded in this way until I had ten informative
genes flanking the original diploid gene, five on each side. I discarded the
syntenic orthogroup if I had to check more than 25 genes on one side, or if I
ran off the chromosome before I had 5 good neighbors.
Next, I checked that at least three of my five neighboring genes on either side 
of the original diploid genes had orthologs in the same 200kb region of the
Bhyb26 genome. 
Finally, I pulled the upstream and downstream neighbors that were closest to
the original diploid gene, and extracted the region between and including these
two neighbor genes. This Bhyb26 In 579 out of 588 cases, all 

Currently I'm telling it to grab the two Bhyb26 orthologs of the 
diploid genes flanking the diploid ortholog of the 'lost' Bhyb26 
gene (wow that was a mouthful). If it can't find a syn ortho in 
Bhyb26 for either (or both) flanking genes, it walks outward along
the chromosome until it does (usually it does not have to go far). 
If the final flanking genes are not within 20kb of each other, give 
up on this 'missing gene'. 

Last run:
664 orthogroups remained after selecting those with a single gene in the appropriate diploid.
588 orthogroups remained after requiring that there be ten flanking (or nearly flanking) diploid genes with a single ortholog in the appropriate Bhyb26 subgenome.
534 orthogroups remained after requiring that at least 4 of the 5 B. hybridum orthologs on either side of the diploid gene be within 200kb of each other.
517 orthogroups remained after requiring that the final upstream and downstream B. hybridum 'anchor' genes demarcating my candidate region be within 200kb of each other.
Ended up with 517 regions representing the Bhyb26 genomic neighborhood where the lost gene should be.
"""

from itertools import compress
from collections import Counter

class Gene(object):
	def __init__(self, genome, geneID, chrom, start, end, strand, synOgID):
		self.genome = genome
		self.geneID = geneID
		self.chrom = chrom
		self.start = int(start)
		self.end = int(end)
		self.strand = strand
		self.synOgID = synOgID

def get_gff_data_structures(gfffile):
	gfflist = []
	gffdic = {}
	with open(gfffile, 'r') as gff:
		counter = 0
		for line in gff.readlines()[3:]:
			L = line.strip('\n').split('\t') 
			if L[2]=='gene':
				mygeneID = L[8].split(';')[1][5:]
				gfflist.append(mygeneID)
				gffdic[mygeneID] = counter
				counter += 1
	return([gfflist, gffdic])

def two_or_more_genes_are_close_together(genelist): #if they aren't, function returns an empty list. 
#if they are, returns a list containing a tuple of the two gene objs in arbitrary order. 
#if more than two genes in genelist are close together, returns all pairs that are close together. 
#returns a non-redundant list, e.g. [ (geneA, geneB) ] not [ (geneA, geneB), (geneB, geneA) ] 
	list_a = [(x,y) for i,x in enumerate(genelist) for j,y in enumerate(genelist) if i > j] #a list of tuples, where each tuple contains a pair of gene objs
	fil = [(x.chrom == y.chrom and abs(int(x.start) - int(y.start)) < kbwindow*1000) for i,x in enumerate(genelist) for j,y in enumerate(genelist) if i > j] #a boolean list
	return(list(compress(list_a, fil)))


def walk_along(dipgeneobj, genesdic, geneslist, has_one): # Searches outward from a starting diploid gene 
#in one direction (first upstream, then downstream) until it finds five diploid genes that each have a 1:1 ortholog 
#If it doesn't find five such genes after searching n genes, gives up
#Note, we could change this to allow e.g. 2:1 relationships, i.e. multiple copies in the diploid but one in the polyploid
	mygeneindex = genesdic[dipgeneobj.geneID]
	upstream = []
	downstream = []
	current_index = mygeneindex
	while len(upstream) < 5:
		current_gene = geneIDs_to_geneobjs[geneslist[current_index]]
		if current_gene.geneID in has_one:
			upstream.append(current_gene)
		current_index -= 1
		if current_index < (mygeneindex - 25) or current_gene.chrom != dipgeneobj.chrom: #give up if we have to move beyond n genes to get 
		#5 genes with a single Bhyb26 ortholog or if we run off the chromosome
			break
	current_index = mygeneindex
	while len(downstream) < 5:
		current_gene = geneIDs_to_geneobjs[geneslist[current_index]]
		if current_gene.geneID in has_one:
			downstream.append(current_gene)
		current_index += 1
		if current_index > (mygeneindex + 25) or current_gene.chrom != dipgeneobj.chrom:
			break
	return(upstream, downstream)

def get_ortholog(dipgeneobj): #expects that you've already checked that the diploid gene has one ortholog
	if dipgeneobj.genome == 'Bdistachyon':
		return(geneIDs_to_geneobjs[has_one_BhD_ortholog[dipgeneobj.geneID]])
	if dipgeneobj.genome == 'Bstacei':
		return(geneIDs_to_geneobjs[has_one_BhS_ortholog[dipgeneobj.geneID]])



################################## START ##################################

#First, create gene objects. Make a dictionary of syntenic orthogroup IDs that point to gene objects,
#and another dictionary of gene IDs that point to syntenic orthogroup IDs. This allows us to easily
#access the genes in a syntenic orthogroup and the gff data for those genes.

kbwindow = 200 #genes are considered "close together" if their start coordinates are less than kbwindow*1000 bp apart. 
synOgdict = {} #key is a synOgID, value is a LIST of Gene objects
geneIDs_to_synOGs = {} #keys are geneIDs, values are synOgIDs
geneIDs_to_geneobjs = {}
with open('gffWithOgs.txt', 'r') as syntenyFile: #made from genespace run GENESPACE_all_genomes
	for line in syntenyFile.readlines()[1:]:
		L = line.strip('\n').split('\t')
		myGene = Gene(L[0], L[1], L[2], L[3], L[4], L[5], L[14])
		geneIDs_to_synOGs[L[1]] = L[14]
		geneIDs_to_geneobjs[L[1]] = myGene
		if L[14] in synOgdict:
			synOgdict[L[14]].append(myGene)
		else:
			synOgdict[L[14]] = [myGene]

#make dictionaries mapping diploid genes to Bhyb26 genes, note these are not syntenic orthologs,
#just orthofinder best hits. I'm keeping one-to-one and many-to-one relationships, just discarding
#many-to-many. It made no difference whether I used orthologs or syntenic orthogroups, tried both ways.
has_one_BhD_ortholog = {} 
has_one_Bd_ortholog = {}
with open('Bdistachyon__v__Bhyb26D.tsv', 'r') as myfile:
	read = myfile.readlines()
	for n in range(1, len(read)):
		line = read[n]
		L = line.strip('\n').split('\t')
		if not "," in L[1]:
			if not "," in L[2]: #one to one
				Bdgene, BhDgene = L[1], L[2]
				has_one_BhD_ortholog[Bdgene] = BhDgene
				has_one_Bd_ortholog[BhDgene] = Bdgene
			if "," in L[2]: #one diploid to many polyploid
				Bdgene, BhDgenes = L[1], L[2].split(", ")
				for BhDgene in BhDgenes:
					has_one_Bd_ortholog[BhDgene] = Bdgene
		if "," in L[1]:
			if not "," in L[2]: #many diploid to one polyploid
				Bdgenes, BhDgene = L[1].split(", "), L[2]
				for Bdgene in Bdgenes:
					has_one_BhD_ortholog[Bdgene] = BhDgene

has_one_BhS_ortholog = {} 
has_one_Bs_ortholog = {}
with open('Bstacei__v__Bhyb26S.tsv', 'r') as myfile:
	read = myfile.readlines()
	for n in range(1, len(read)):
		line = read[n]
		L = line.strip('\n').split('\t')
		if not "," in L[1]:
			if not "," in L[2]: #one to one
				Bsgene, BhSgene = L[1], L[2]
				has_one_BhS_ortholog[Bsgene] = BhSgene
				has_one_Bs_ortholog[BhSgene] = Bsgene
			if "," in L[2]: #one diploid to many polyploid
				Bsgene, BhSgenes = L[1], L[2].split(", ")
				for BhSgene in BhSgenes:
					has_one_Bs_ortholog[BhSgene] = Bsgene
		if "," in L[1]:
			if not "," in L[2]: #many diploid to one polyploid
				Bsgenes, BhSgene = L[1].split(", "), L[2]
				for Bsgene in Bsgenes:
					has_one_BhS_ortholog[Bsgene] = BhSgene

#Next, make a list of genes ordered by their order on the chromosome, and make a dictionary of
#gene IDs pointing to their index in the ordered list. Do this for each diploid progenitor.
#We will use these to find the chromosomal neighborhood of the 'missing' Bhyb26 genes.

Bdgffstuff = get_gff_data_structures('Bdistachyon_556_v3.2.gene.gff3')
Bdgenesordered = Bdgffstuff[0] #a list of genes ordered by their order on the chromosome
Bdgenesdic = Bdgffstuff[1] #keys are geneIDs, values are their order in Bdgenesordered, i.e. their order in the genome

Bsgffstuff = get_gff_data_structures('Bstacei_316_v1.1.gene.gff3')
Bsgenesordered = Bsgffstuff[0]
Bsgenesdic = Bsgffstuff[1] #keys are geneIDs, values are their order in Bsgenesordered, i.e. their order in the genome

Bhgffstuff = get_gff_data_structures('BhybridumBhyb26v2.1.gene.gff3')
Bhgenesordered = Bhgffstuff[0]
Bhgenesdic = Bhgffstuff[1] #keys are geneIDs, values are their order in Bhgenesordered, i.e. their order in the genome


#Get orthogroups with 'missing' genes

mysynOgs = {} #keys are orthogroupIDs, value is the subgenome from which a gene has been putatively lost
with open('lostgenes.txt', 'r') as synOgfile: #made with my script bin_orthogroups.R
	for line in synOgfile.readlines()[1:]:
		L = line.strip('\n').split('\t')
		mysynOgs[L[0]] = L[1]

#Iterate through the orthogroups of interest. Keep only those that have ONE diploid
#gene corresponding to the missing subgenome. For example, if the missing gene is
#from BhD and the orthogroup has two Bd genes, we don't know which one to use as the
#best ortholog. It's possible that those two Bd genes are tandem duplicates and are
#syntenic to the same region of BhD, but I'm not going to try to rescue possible 
#tandem duplicates at this point, it's too much effort for too little reward.

#Then get the ten diploid genes flanking the diploid ortholog: five on each side.
#Next, we run the same test on each of the flanking genes: does it have a single
#syntenic ortholog in the appropriate subgenome? If not, keep 'walking' up or down
#the chromosome until we have five usable genes on each side of the original diploid ortholog.

rnd0_pass = [] #synOgIDs for synOgs that had one and only one gene in the appropriate diploid
rnd1_pass = [] #synOgIDs for synOgs for which I found ten flanking genes
upstream_flanking_genes = {} #key is synOg ID, value is a list of gene objects
downstream_flanking_genes = {} 
orig_dip_genes = {} #key is synOg ID, value is a gene object for the central diploid gene
for mysynOg, subg in mysynOgs.items():
	upstream_dipgenes = []
	downstream_dipgenes = []
	if subg == 'Bhyb26D':
		if len([g for g in synOgdict[mysynOg] if g.genome=='Bdistachyon']) == 1:
			rnd0_pass.append(mysynOg)
			original_dipgene = [g for g in synOgdict[mysynOg] if g.genome=='Bdistachyon'][0]
			#walk_along(dipgeneobj, genesdic, geneslist, has_one)
			upstream_dipgenes, downstream_dipgenes = walk_along(original_dipgene, Bdgenesdic, Bdgenesordered, has_one_BhD_ortholog) 
	elif subg == 'Bhyb26S':
		if len([g for g in synOgdict[mysynOg] if g.genome=='Bstacei']) == 1:
			rnd0_pass.append(mysynOg)
			original_dipgene = [g for g in synOgdict[mysynOg] if g.genome=='Bstacei'][0]
			#walk_along(dipgeneobj, genesdic, geneslist, has_one)
			upstream_dipgenes, downstream_dipgenes = walk_along(original_dipgene, Bsgenesdic, Bsgenesordered, has_one_BhS_ortholog) 
	if len(upstream_dipgenes) == 5 and len(downstream_dipgenes) == 5: # we may not get ten flanking genes, e.g. if synOg is in/near a telomere
		rnd1_pass.append(mysynOg)
		upstream_flanking_genes[mysynOg] = upstream_dipgenes 
		downstream_flanking_genes[mysynOg] = downstream_dipgenes
		orig_dip_genes[mysynOg] = original_dipgene

#Now that we have ten flanking genes, we need to find the candidate pseudogene-containing
#region in B. hybridum. Average gene length in B. hybridum is about 3.5 kb.
#Start by comparing each of the flanking genes to all the others. Basically, we will
#identify all the Bhyb26 "neighborhoods" that the ten genes point to. If two genes'
#Bhyb26 orthologs are within 200kb of each other, we will say they belong to the same 
#neighborhood. (I'm deliberately being generous). Then we will choose the neighborhood 
#with the highest representation among our flanking genes, and discard any diploid genes 
#that don't point to that neighborhood. Finally, from those flanking genes that remain, we 
#will choose the two most central ones to be our "anchors", and the region between them 
#will be the candidate pseudogene region.


rnd2_pass = [] #synOgIDs for synOgs for which I found at least four genes on either side pointed to the same Bhyb26 region (as each other)
final_upstream_flanking_genes = {} #key is synOg ID, value is a list of gene objects
final_downstream_flanking_genes = {} 
for mysynOg in rnd1_pass:
	upstream_orthos = [get_ortholog(g) for g in upstream_flanking_genes[mysynOg]] 
	downstream_orthos = [get_ortholog(g) for g in downstream_flanking_genes[mysynOg]]
	upstream_all_neighbors = two_or_more_genes_are_close_together(list(set(upstream_orthos))) #unique-ify in case multiple neighbors point to the same ortholog
	downstream_all_neighbors = two_or_more_genes_are_close_together(list(set(downstream_orthos)))
	up_flat_list = [gene for sublist in upstream_all_neighbors for gene in sublist]
	down_flat_list = [gene for sublist in downstream_all_neighbors for gene in sublist]
	up_flat_list = [tup[0] for tup in Counter(up_flat_list).most_common()] 
	down_flat_list = [tup[0] for tup in Counter(down_flat_list).most_common()] 
	#get the maximum number of times any gene appears, and all genes that appear that number of times.
	if len(up_flat_list) >= 4 and len(down_flat_list) >= 4: #this discards "neighborhoods" with less support
		final_upstream_flanking_genes[mysynOg] = up_flat_list
		final_downstream_flanking_genes[mysynOg] = down_flat_list
		rnd2_pass.append(mysynOg)

#BUG TESTING
# I THINK THERE IS AN ISSUE WITH two_or_more_genes_are_close_together()
# upstream_orthos = [get_ortholog(g) for g in upstream_flanking_genes[mysynOg]]
# upstream_all_neighbors = two_or_more_genes_are_close_together(upstream_orthos)
# up_flat_list = [gene for sublist in upstream_all_neighbors for gene in sublist]
#  = [tup[0] for tup in Counter(up_flat_list).most_common()] 
# downstream_orthos = [get_ortholog(g) for g in downstream_flanking_genes[mysynOg]] 
# downstream_all_neighbors = two_or_more_genes_are_close_together(list(set(downstream_orthos))) #unique-ify in case multiple neighbors point to the same ortholog
# down_flat_list = []
# for tup in downstream_all_neighbors:
# 	down_flat_list.append(tup[0])
# 	down_flat_list.append(tup[1])
# [gene for sublist in downstream_all_neighbors for gene in sublist] 
#  = [tup[0] for tup in Counter(down_flat_list).most_common()] 
# #get the maximum number of times any gene appears, and all genes that appear that number of times.
# #this discards "neighborhoods" with less support
# if len(up_flat_list) >= 4 and len(down_flat_list) >= 4:
# 	final_upstream_flanking_genes[mysynOg] = up_flat_list
# 	final_downstream_flanking_genes[mysynOg] = down_flat_list
# 	rnd2_pass.append(mysynOg)

# bugtesting: 
# cnt=Counter()
# myL = []
# for og in final_upstream_flanking_genes.values():
#  	myL.append(len(og))

# for mylen in myL:
# 	cnt[mylen]+=1
#  
# cnt
# Counter({5: 458, 4: 102, 3: 19, 2: 6})
# Downstream info: Counter({5: 449, 4: 112, 3: 15, 2: 9})

# Convert our diploid flanking genes to their hybridum orthologs.
# Get the upstream HYBRIDUM gene with the highest start coord and the downstream HYBRIDUM gene with the lowest stop coord.
# Our final candidate pseudogene region will be the region between and including these two genes. We are getting the
# two most central genes we can so that our candidate region will be as small as possible. 

final_flank = {} #key is a synOgID, value is a list: [upstream gene, downstream gene] (Gene objects, that is, from HYBRIDUM)	
for mysynOg in rnd2_pass:
	highest_coord = max([g.start for g in final_upstream_flanking_genes[mysynOg]]) 
	for upgene in final_upstream_flanking_genes[mysynOg]:
		if upgene.start == highest_coord or upgene.end == highest_coord:
			final_flank[mysynOg] = [upgene] # get the ortholog of the upstream diploid gene that is as close as possible to the original diploid gene
			break 
	lowest_coord = min([g.end for g in final_downstream_flanking_genes[mysynOg]]) 
	for downgene in final_downstream_flanking_genes[mysynOg]:
		if downgene.start == lowest_coord or downgene.end == lowest_coord:
			final_flank[mysynOg].append(downgene) # get the ortholog of the downstream diploid gene that is as close as possible to the original diploid gene
			break 

#Finally, get the gff info for the two genes that mark the borders of our region putatively containing the candidate pseudogene.
#There will be cases where the syn orthos are in the reverse order in B. hybridum, i.e. the syn ortho of the downstream
#diploid gene is upstream of the syn ortho of the upstream diploid gene, presumably due to an inversion 
#(I quickly scanned these cases by eye and the diploid genes are very close together in all cases)
#in these cases, still grab the area between hyb1 and hyb2, just go from hyb2.start to hyb1.end. 
#(Drew this out with pen and paper and I'm quite sure it works out)

hybridum_chroms = []
hybridum_starts = []
hybridum_ends = []
diploid_genes = []
rnd3_pass = []
for mysynOg in rnd2_pass:
	hyb1 = final_flank[mysynOg][0]
	hyb2 = final_flank[mysynOg][1]
	gene_interval_start = hyb1.start
	gene_interval_end = hyb2.end
	if gene_interval_start > gene_interval_end:
		gene_interval_start = hyb2.start
		gene_interval_end = hyb1.end
	if gene_interval_end - gene_interval_start < kbwindow*1000: #The most egregious cases that don't pass here (region is > 1Mb) are all cases 
	#where a structural rearrangement in Bhyb26 placed the upstream and downstream orthologs far apart from each other. (I checked)
	#which is kind of interesting, because it suggests the missing gene is smack on the border of a rearrangement. There were 5 such cases.
		rnd3_pass.append(mysynOg)
		diploid_genes.append(orig_dip_genes[mysynOg].geneID)
		hybridum_chroms.append(hyb1.chrom)
		hybridum_starts.append(str(gene_interval_start))
		hybridum_ends.append(str(gene_interval_end))

	# else: #bug testing
	# 	print("Final upstream hybridum gene coordinates:")
	# 	print([g.start for g in final_upstream_flanking_genes[mysynOg]])
	# 	print("Final downstream hybridum gene coordinates:")
	# 	print([g.start for g in final_downstream_flanking_genes[mysynOg]])

		#print("NOTE! Best upstream and downstream genes are quite far apart for synOg {}!".format(mysynOg))
		#print("The entire candidate pseudogene region will be {}kb.".format(
		#	round((gene_interval_end - gene_interval_start)/1000) ))




print("Started with {} orthogroups containing a putative lost gene.".format(len(mysynOgs)))
print("{} orthogroups remained after selecting those with a single gene in the appropriate diploid.".format(len(rnd0_pass)))
print("{} orthogroups remained after requiring that there be ten flanking (or nearly flanking) diploid genes with a single ortholog in the appropriate Bhyb26 subgenome.".format(len(rnd1_pass)))
print("{} orthogroups remained after requiring that at least 4 of the 5 B. hybridum orthologs on either side of the diploid gene be within {}kb of each other.".format(len(rnd2_pass), kbwindow))
print("{} orthogroups remained after requiring that the final upstream and downstream B. hybridum 'anchor' genes demarcating my candidate region be within {}kb of each other.".format(len(rnd3_pass), kbwindow))
print("Ended up with {} regions representing the Bhyb26 genomic neighborhood where the lost gene should be.".format(len(hybridum_starts)))



with open('missing_gene_intervals.bed', 'w') as outF:
	for n in range(len(diploid_genes)):
		outF.write('{}\t{}\t{}\t{}\n'.format(hybridum_chroms[n], hybridum_starts[n], hybridum_ends[n], diploid_genes[n]))





