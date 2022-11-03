import glob
import os
import shutil

def parse_conserved_gene_info():
	mydata = {}
	with open('conservedgenes.txt', 'r') as myfile:
		read = myfile.readlines()
		for n in range(1, len(read)):
			L = read[n].split('\t')
			if len(L) == 4: #two progenitors
				PP1, PP2, dipgene1, dipgene2 = L[0], L[1], L[2], L[3].strip('\n')
				mydata[PP1] = dipgene1
				mydata[PP2] = dipgene2
			elif len(L) == 3: #one progenitor
				PP1, PP2, dipgene = L[0], L[1], L[2].strip('\n')
				mydata[PP1] = dipgene
				mydata[PP2] = dipgene
	return(mydata)


PPtodiploid = {}
if os.path.exists('conservedgenes.txt'):
	PPtodiploid = parse_conserved_gene_info()

fafilelist = glob.glob('candidatePseudogeneRegions/*.region.fa')
#If we are working on missing genes, these will be diploid gene IDs.
#If we are working on conserved genes, these will be polyploid gene IDs.
#I did this because the gene IDs must be unique.

#Check whether we have old files left over from previous runs of this script.
#We will delete these and get a fresh start.
for f in fafilelist:
	if f[27:-10].isdigit():
		os.remove(f)

fafilelist = [f for f in fafilelist if not f[27:-10].isdigit()]
# f[27:-10] strips away the 'candidatePseudogeneRegions/' and '.region.fa'
# If the fasta file is e.g. 33.region.fa, that means we have run this script previously.
		
if PPtodiploid: #if we are working with conserved genes:
	with open('aliases.txt', 'w') as outF:
		for n in range(len(fafilelist)):
			PPgeneID = fafilelist[n][27:-10] #strip away the 'candidatePseudogeneRegions/' and '.region.fa'
			alias = str(n+1)
			outF.write('{}\t{}\n'.format(alias, PPgeneID))
			shutil.copy('candidatePseudogeneRegions/{}.region.fa'.format(PPgeneID), 'candidatePseudogeneRegions/{}.region.fa'.format(alias))
			#There are two polyploid genes for each diploid gene, so we will copy and rename the peptide file twice.
			#This means the same file will appear in two different names, e.g. 145.pep.fa and 2067.pep.fa might be identical,
			#so we will have a lot of files, but I think it will be nicely organized to have 1 alias : 1 region file : 1 peptide file.
			shutil.copy('candidatePseudogeneRegions/{}.pep.fa'.format(PPtodiploid[PPgeneID]), 'candidatePseudogeneRegions/{}.pep.fa'.format(alias))
			os.remove('candidatePseudogeneRegions/{}.region.fa'.format(PPgeneID)) #make sure we've successfully created the new files before deleting the old ones
	for dipgene in set(PPtodiploid.values()): #clean up the peptide files we no longer need (since we've copied and renamed them)
		os.remove('candidatePseudogeneRegions/{}.pep.fa'.format(dipgene))
		

if not PPtodiploid: #if we are working with missing genes:
	with open('aliases.txt', 'w') as outF:
		for n in range(len(fafilelist)):
			dipgeneID = fafilelist[n][27:-10] #strip away the 'candidatePseudogeneRegions/' and '.region.fa'
			alias = str(n+1)
			outF.write('{}\t{}\n'.format(str(n+1), dipgeneID))
			#There is only one polyploid gene per diploid gene, so we can immediately rename the peptide files.
			os.rename('candidatePseudogeneRegions/{}.region.fa'.format(dipgeneID), 'candidatePseudogeneRegions/{}.region.fa'.format(alias))
			os.rename('candidatePseudogeneRegions/{}.pep.fa'.format(dipgeneID), 'candidatePseudogeneRegions/{}.pep.fa'.format(alias))