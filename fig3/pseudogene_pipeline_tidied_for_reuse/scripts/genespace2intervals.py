"""
A script to follow up on my 'missing' genes.

This script takes as input gffWithOgs.txt from GENESPACE and
a file 'lost_genes.txt', produced by my script bin_orthogroups.R.

The output is a bed file containing the coordinates of genomic regions
that should contain a pseudogene. The diploid ortholog is listed in the 'names' column.

The workflow is as follows: 
First, syntenic orthogroups were filtered to discard any that orthogroups
that had more than one gene from the diploid corresponding to the missing 
subgenome. (In my experiment, none of the orthogroups had this issue.) 
Once I had identified the diploid gene corresponding to the 'missing' Bhyb26 
gene, I 'walked' outward along the diploid chromosome in both directions, checking 
whether each nearby gene had a single ortholog in the appropriate polyploid 
subgenome. If it had no orthologs or many orthologs, I skipped it and moved 
on to the next-closest gene. I proceeded in this way until I had ten informative
genes flanking the original diploid gene, five on each side. I discarded the
syntenic orthogroup if I had to check more than 25 genes on one side, or if I
reached the end of the chromosome before I had 5 good neighbors.
Next, I checked that at least four of my five neighboring genes on either side 
of the original diploid genes had orthologs in the same 200kb region of the
polyploid genome. In other words, 8 of my 10 valid diploid neighbors had to point
to the same genomic region in the polyploid. 
Finally, I pulled the upstream and downstream neighbors that were closest to
the original diploid gene, and extracted the region between and including these
two neighbor genes. 
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
from collections import Counter, defaultdict

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

def walk_along(dipgeneobj, genesdic, geneslist): # Searches outward from a starting diploid gene 
#in one direction (first upstream, then downstream) until it finds five diploid genes that have a 1:1 ortholog 
#If it doesn't find five such genes after searching <max_gene_search> genes, gives up
	mygeneindex = genesdic[dipgeneobj.geneID]
	upstream = []
	downstream = []
	current_index = mygeneindex
	while len(upstream) < 5:
		current_gene = geneIDs_to_geneobjs[geneslist[current_index]]
		if current_gene.geneID in orthology_dict:
			if len(orthology_dict[current_gene.geneID]) == 1:
				upstream.append(current_gene)
		current_index -= 1
		if current_index < (mygeneindex - max_gene_search) or current_gene.chrom != dipgeneobj.chrom: #give up if we have to move 
			#beyond n genes to get 5 genes with a single polyploid ortholog, or if we run off the chromosome
			break
	current_index = mygeneindex
	while len(downstream) < 5:
		current_gene = geneIDs_to_geneobjs[geneslist[current_index]]
		if current_gene.geneID in orthology_dict:
			if len(orthology_dict[current_gene.geneID]) == 1:
				downstream.append(current_gene)
		current_index += 1
		if current_index > (mygeneindex + max_gene_search) or current_gene.chrom != dipgeneobj.chrom:
			break
	return(upstream, downstream)

def get_ortholog(dipgeneobj): #expects that you've already checked that the diploid gene has one ortholog
	if dipgeneobj.genome == 'Bdistachyon':
		return(geneIDs_to_geneobjs[has_one_BhD_ortholog[dipgeneobj.geneID]])
	if dipgeneobj.genome == 'Bstacei':
		return(geneIDs_to_geneobjs[has_one_BhS_ortholog[dipgeneobj.geneID]])

def two_or_more_genes_are_close_together(genelist): #if they aren't, function returns an empty list. 
#if they are, returns a list containing a tuple of each pair of 'close together' gene objs in arbitrary order. 
#if more than two genes in genelist are close together, returns all pairs that are close together. 
#returns a non-redundant list, e.g. [ (geneA, geneB) ] not [ (geneA, geneB), (geneB, geneA) ] 
	list_a = [(x,y) for i,x in enumerate(genelist) for j,y in enumerate(genelist) if i > j] #a list of tuples, where each tuple contains a pair of gene objs
	fil = [(x.chrom == y.chrom and abs(int(x.start) - int(y.start)) < kbwindow*1000) for i,x in enumerate(genelist) for j,y in enumerate(genelist) if i > j] #a boolean list
	return(list(compress(list_a, fil)))

#example input for two_or_more_genes_are_close_together():
#[<__main__.Gene object at 0x7f8d71d3bac0>, <__main__.Gene object at 0x7f8d71d3ba90>, <__main__.Gene object at 0x7f8d71d3bbe0>, <__main__.Gene object at 0x7f8d71d3ba90>]
#The corresponding geneIDs for that would be:
#['BrBhyv21016157m.g', 'BrBhyv21016147m.g', 'BrBhyv21016146m.g', 'BrBhyv21016145m.g'] 
#^Only 4 genes because two of the valid diploid genes pointed to the same polyploid ortholog

#example result from two_or_more_genes_are_close_together():
#[ (<__main__.Gene object at 0x7f8d71d3bac0>, <__main__.Gene object at 0x7f8d71d3ba90>), 
#	(<__main__.Gene object at 0x7f8d71d3bbe0>, <__main__.Gene object at 0x7f8d71d3ba90>), 
#	(<__main__.Gene object at 0x7f8d71d3bbe0>, <__main__.Gene object at 0x7f8d71d3bac0>), 
#	(<__main__.Gene object at 0x7f8d71d3baf0>, <__main__.Gene object at 0x7f8d71d3ba90>), 
#	(<__main__.Gene object at 0x7f8d71d3baf0>, <__main__.Gene object at 0x7f8d71d3bac0>), 
#	(<__main__.Gene object at 0x7f8d71d3baf0>, <__main__.Gene object at 0x7f8d71d3bbe0>) ]

#The corresponding geneIDs for that would be:
# [ ('BrBhyv21016146m.g', 'BrBhyv21016145m.g'),
# ('BrBhyv21016157m.g', 'BrBhyv21016145m.g'),
# ('BrBhyv21016157m.g', 'BrBhyv21016146m.g'),
# ('BrBhyv21016147m.g', 'BrBhyv21016145m.g'),
# ('BrBhyv21016147m.g', 'BrBhyv21016146m.g'),
# ('BrBhyv21016147m.g', 'BrBhyv21016157m.g') ]

#In most cases, there will be many pairs of genes, indicating that several genes in the input list are close together.
#In this example, all of the polyploid genes were close together.


################################## START ##################################


######################### EDIT HERE #########################
kbwindow = 200 #genes are considered "close together" if their start coordinates are less than kbwindow*1000 bp apart. 
max_gene_search = 25 #give up searching for diploid genes with a single polyploid ortholog after we have searched n genes in one direction
polyploids = ['Bhyb26D', 'Bhyb26S']
progenitors = ['Bdistachyon', 'Bstacei'] #if you have only one progenitor, just list it TWICE
#^These MUST be in corresponding order! polyploids[0] should correspond to progenitors[0]
progenitor_gffs = ['Bdistachyon_556_v3.2.gene.gff3', 'Bstacei_316_v1.1.gene.gff3']
#^Order matters here, too. Again, if only one progenitor, just list the same file twice 
path2gffwithOGs = '/home/virginia/Documents/school/vogelLab/GENESPACE_test/results/gffWithOgs.txt' #from GENESPACE
path2lostgenes = '/home/virginia/Documents/school/vogelLab/notebook/2022/lost_genes/lost_genes/lostgenes.txt' #made with my script bin_orthogroups.R
ortholog_files = ['Bdistachyon__v__Bhyb26D.tsv', 'Bstacei__v__Bhyb26S.tsv'] #from orthofinder via GENESPACE
#^^^Need one ortholog file for each polyploid subgenome. Diploid__v__Polyploid.tsv
# With any of these files, you can either include the full path to the file elsewhere
# or just keep the file in your working directory--it shouldn't matter.
#############################################################


#First, create som useful data structures. 
geneIDs_to_geneobjs = {} #key = the geneID (a string), value = the corresponding instance of class Gene
synOgdict = {} #key is a synOg ID, value is a dictionary. Within that dict, key is the genome (a string), value is a list containing gene ID(s).
#e.g. { 
#	'21902':{'ABR113D':['Brahy.D01G0000100'], 'ABR113S':['Brahy.S02G0423800'], 'Bstacei':['Brast02G394300'], 'Bdistachyon':['Bradi1g00485']},
#	'68269':{'Osativa':['LOC_Os06g46690'], 'Bstacei':['Brast02G106200', 'Brast02G106300']},
#	...
#	}
with open(path2gffwithOGs, 'r') as syntenyFile: 
	for line in syntenyFile.readlines()[1:]:
		L = line.strip('\n').split('\t')
		genome, geneID, synOg = L[0], L[1], L[13] #IMPORTANT! BE CAREFUL!!! Make sure the column number for synOg is correct (zero-based counting).
		#^The columns have changed with different versions of GENESPACE.
		myGene = Gene(genome, geneID, L[2], L[3], L[4], L[5], synOg) 
		geneIDs_to_geneobjs[L[1]] = myGene
		if synOg not in synOgdict:
			synOgdict[synOg] = defaultdict(list) #default dict simply allows us to append to a list even if it doesn't exist yet. (The list will be created.)
		synOgdict[synOg][genome].append(geneID)


#Get orthogroups with 'missing' genes
mysynOgs = {} #keys are synOgIDs, value is the subgenome from which a gene has been putatively lost
with open(path2lostgenes, 'r') as synOgfile: 
	for line in synOgfile.readlines()[1:]:
		L = line.strip('\n').split('\t')
		mysynOgs[L[0]] = L[1]


orthology_dict = {} # I believe these are a bit more accurate than the syntenic orthologs.
# ^ Here, keys are diploid and polyploid gene IDs, and value is a list of gene objects:
# that gene's orthologs in the corresponding diploid or polyploid subgenome.
for f in ortholog_files:
	with open(f, 'r') as orthofile:
		read = orthofile.readlines()
		for n in range(1, len(read)):
			L = read[n].split('\t')
			dipgenes, PPgenes = L[1].split(', '), L[2].strip('\n').split(', ') #lists with at least one element
			for g in dipgenes:
				orthology_dict[g] = [geneIDs_to_geneobjs[gID] for gID in PPgenes]
			for g in PPgenes:
				orthology_dict[g] = [geneIDs_to_geneobjs[gID] for gID in dipgenes]



#Next, get two more data structures: 
#progenitor_genes_ordered is a list of genes ordered by their order on the chromosome.
#progenitor_gene_indices is a dictionary, In that dictionary, keys are geneIDs, and values are their index in the list.
#For Bmex, with only one progenitor, the info for progenitor 1 and progenitor 2 will be the same.

progenitor1_genes_ordered, progenitor1_gene_indices = get_gff_data_structures(progenitor_gffs[0])
progenitor2_genes_ordered, progenitor2_gene_indices = get_gff_data_structures(progenitor_gffs[1])


#Iterate through the diploid genes. Keep only those that have ONE diploid
#gene corresponding to the missing subgenome. For example, if the missing gene is
#from BhD and the orthogroup has two Bd genes, we don't know which one to use as the
#best ortholog. 

#Then get the ten diploid genes flanking the diploid ortholog: five on each side.
#Each of the flanking genes must also have a single syntenic ortholog in the appropriate subgenome.
#Keep 'walking' up or down the chromosome until we have five valid genes on each side of the original diploid ortholog.


rnd0_pass = [] #synOgs that had one and only one gene in the appropriate diploid
rnd1_pass = [] #synOgs for which I found ten valid flanking genes
upstream_flanking_genes = {} #key is synOg ID, value is a list of gene objects
downstream_flanking_genes = {} 
orig_dip_genes = {} #key is synOg ID, value is a gene object for the central diploid gene
for mysynOg, subg in mysynOgs.items():
	myprog = progenitors[polyploids.index(subg)] # name of progenitor that corresponds to the missing subgenome, e.g. 'Bdistachyon'
	upstream_dipgenes, downstream_dipgenes = [], []
	if len(synOgdict[mysynOg][myprog]) == 1: #Does this synOg have only one gene from this diploid progenitor? Should be true in most cases.
		rnd0_pass.append(mysynOg)
		mygeneobj = geneIDs_to_geneobjs[ synOgdict[mysynOg][myprog][0] ]
		if myprog == progenitors[0]: #if progenitor1 == progenitor2, then this 'if' part of the 'if else' statement will always be executed
			upstream_dipgenes, downstream_dipgenes = walk_along(mygeneobj, progenitor1_gene_indices, progenitor1_genes_ordered)
		else:
			upstream_dipgenes, downstream_dipgenes = walk_along(mygeneobj, progenitor2_gene_indices, progenitor2_genes_ordered)
		if len(upstream_dipgenes) == 5 and len(downstream_dipgenes) == 5: # Did we manage to get ten valid flanking genes?
			rnd1_pass.append(mysynOg)
			upstream_flanking_genes[mysynOg] = upstream_dipgenes 
			downstream_flanking_genes[mysynOg] = downstream_dipgenes
			orig_dip_genes[mysynOg] = mygeneobj




#Now that we have ten polyploid genes whose progenitor orthologs are in the same genomic region, 
#we need to check whether those polyploid genes are also close together. If they are, 
#we will have our candidate pseudogene-containing region in the polyploid. 

#Start by comparing each of the flanking genes to all the others. If two genes are within 200kb 
#of each other in the polyploid genome, we will say they belong to the same 
#neighborhood. (I'm deliberately being generous). 

#two_or_more_genes_are_close_together essentially gives us a list of neighborhoods.
#For example, the result
#[(geneA, geneB), (geneA, geneB, geneC)] tells us that there are two 200kb neighborhoods
#that contain at least two of our polyploid genes, and one of those neighborhoods contains three.

#Then we will choose a neighborhood 
#with the highest representation among our flanking genes, and discard any diploid genes 
#that don't point to that neighborhood. Finally, from those flanking genes that remain, we 
#will choose the two most central ones to be our "anchors", and the region between them 
#will be the candidate pseudogene region.


rnd2_pass = [] #synOgs for which I found at least four genes on either side pointing to the same polyploid region 
final_upstream_flanking_genes = {} #key is synOg ID, value is a list of gene objects
final_downstream_flanking_genes = {} 
for mysynOg in rnd1_pass:
	upstream_orthos = [ orthology_dict[g.geneID][0] for g in upstream_flanking_genes[mysynOg] ]
	#^Now we are getting the polyploid genes (a list of Gene objects). Safe to get the first ([0]) ortholog because we already made sure there is only one.
	downstream_orthos = [ orthology_dict[g.geneID][0] for g in downstream_flanking_genes[mysynOg] ]
	upstream_all_neighbors = two_or_more_genes_are_close_together(list(set(upstream_orthos))) #unique-ify in case multiple neighbors point to the same ortholog
	downstream_all_neighbors = two_or_more_genes_are_close_together(list(set(downstream_orthos)))
	up_flat_list = [gene for sublist in upstream_all_neighbors for gene in sublist] #convert the nested list to a flat list
	down_flat_list = [gene for sublist in downstream_all_neighbors for gene in sublist]
	#Next, get the maximum number of times any gene appears, and keep only the genes that appear that number of times.
	#If all polyploid genes were close together, then this list will be the same as the input for two_or_more_genes_are_close_together().
	#(This will usually be the case.)
	up_flat_list = [tup[0] for tup in Counter(up_flat_list).most_common()] #a list of gene objects. 
	down_flat_list = [tup[0] for tup in Counter(down_flat_list).most_common()] 
	if len(up_flat_list) >= 4 and len(down_flat_list) >= 4: #We will only keep synOgs where at least 4 of 5 flanking genes are close together.
		#NOTE TO SELF: I'm not sure if we are losing some genes at this step if more than 2 diploid genes have the same polyploid ortholog?
		final_upstream_flanking_genes[mysynOg] = up_flat_list
		final_downstream_flanking_genes[mysynOg] = down_flat_list
		rnd2_pass.append(mysynOg)



# Get the upstream polyploid gene with the highest start coord and the downstream polyploid gene with the lowest stop coord.
# Our final candidate pseudogene region will be the region between and including these two genes. 
final_flank = {} #key is a synOgID, value is a list: [upstream gene, downstream gene] (Gene objects)	
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
#There will be cases where the syn orthos are in the reverse order in the polyploid, i.e. the syn ortho of the downstream
#diploid gene is upstream of the syn ortho of the upstream diploid gene, presumably due to an inversion 
#(I quickly scanned these cases by eye and the diploid genes are very close together in all cases)
#in these cases, still grab the area between pp_gene1 and pp_gene2 (polyploid gene), just go from pp_gene2.start to pp_gene1.end. 
polyploid_chroms = []
polyploid_starts = []
polyploid_ends = []
diploid_genes = []
rnd3_pass = []
for mysynOg in rnd2_pass:
	pp_gene1 = final_flank[mysynOg][0]
	pp_gene2 = final_flank[mysynOg][1]
	gene_interval_start = pp_gene1.start
	gene_interval_end = pp_gene2.end
	if gene_interval_start > gene_interval_end:
		gene_interval_start = pp_gene2.start
		gene_interval_end = pp_gene1.end
	if gene_interval_end - gene_interval_start < kbwindow*1000: #The most egregious cases that don't pass here are all cases 
	#where a structural rearrangement in Bhyb26 placed the upstream and downstream orthologs far apart from each other. (I checked)
	#which is kind of interesting, because it suggests the missing gene is smack on the border of a rearrangement. There were 5 such cases.
		rnd3_pass.append(mysynOg)
		diploid_genes.append(orig_dip_genes[mysynOg].geneID)
		polyploid_chroms.append(pp_gene1.chrom)
		polyploid_starts.append(str(gene_interval_start))
		polyploid_ends.append(str(gene_interval_end))


print("Run completed successfully.")
print("Started with {} orthogroups containing a putative lost gene.".format(len(mysynOgs)))
print("{} orthogroups remained after selecting those with a single gene in the appropriate diploid.".format(len(rnd0_pass)))
print("{} orthogroups remained after requiring that there be ten flanking (or nearly flanking) diploid genes with a single ortholog in the appropriate polyploid subgenome.".format(len(rnd1_pass)))
print("{} orthogroups remained after requiring that at least 4 of the 5 B. hybridum orthologs on either side of the diploid gene be within {}kb of each other.".format(len(rnd2_pass), kbwindow))
print("{} orthogroups remained after requiring that the final upstream and downstream 'anchor' genes demarcating my candidate region be within {}kb of each other.".format(len(rnd3_pass), kbwindow))
print("Ended up with {} regions representing the genomic neighborhood where the lost gene should be.".format(len(polyploid_starts)))



with open('missing_gene_intervals.bed', 'w') as outF:
	for n in range(len(diploid_genes)):
		outF.write('{}\t{}\t{}\t{}\n'.format(polyploid_chroms[n], polyploid_starts[n], polyploid_ends[n], diploid_genes[n]))