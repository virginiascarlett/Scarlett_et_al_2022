This is a guide to the latest version of my pseudogene identification protocol.
All my scripts are on my PC in ~/Documents/school/vogelLab/scripts/pseudogenes_redux
(They're on github, too.) This protocol is quite different from (better than) the 2018 version.
Data files are in cori htar archive.

The basic idea is that I took John Lovell's orthogroups and identified those where all
but one (sub)genome was represented. In other words, there was at least one gene from rice, 
both ABR113 subgenomes, and both diploids, but only one of the two Bhyb26 subgenomes. 
I am calling these my 'lost' or 'missing' Bhyb26 genes. I want to investigate those 
chromosomal regions more closely to see whether the gene simply wasn't annotated,
but a fragmented ortholog is still there. If it is still there, we can investigate 
the sequence for signs of pseudogenization/degradation.  

(On my PC, the GENESPACE analysis is in /home/virginia/Documents/school/vogelLab/GENESPACE_all_genomes/.
It's too big for github.) 

It took me a long time to figure out which pairwise comparisons we want to use to
get dN/dS values. For a while, I was thinking of comparing all my brachy genes to 
rice, but that didn't pan out--any differences between the conserved genes and the 
pseudogenes are just swamped by the brachy-rice divergence.

I also spent a lot of time trying to estimate the "time of inactivation", or the time
when a given pseudogene became inactive, based on the method of Chou et al. 2002 
"Inactivation of CMP-N-acetylneuraminic acid hydroxylase occurred prior to brain expansion during human evolution"
Though better explained by Stedman et al. 2004 "Myosin gene mutation correlates with anatomical changes in the human lineage".
A variation of this was implemented in the Xenopus laevis genome paper Session et al. 2016.
I tried plugging in a bunch of different combinations of pairwise dN/dS values into the
formula of Chou et al. 2012. I tried using rice as a reference or using the progenitor
orthologs, I even tried using the homeologs to get my "constrained evolution" dN/dS and
using a universal substitution rate as the dS in their formula. What looked the best
by far was to compare each gene to its diploid ortholog to get the dN/dS, and then use
the conserved genes vs. their orthologs for the 'constrained' dN/dS and the lost genes
vs. their orthologs for the dN of interest. However, even this still looked pretty funky
(see the plot from Mar. 23 2022 in my notebook), even though I'm pretty confident in my
alignments and in my implementation of their formula (see dNdSsignificance.R). We just got
a lot of genes with a negative pseudogenization time, indicating that the dN/dS for the 
pseudogene was actually smaller (more constrained) than for the average conserved gene. 
I ended up concluding that this methodology just isn't a good fit for our system. 

The plot that ended up in the paper (hopefully, as it's still just a manuscript draft as I
write this) was produced by my script dNdSsignificance_all_conserved_genes.R.
How I obtained metrics on premature stop codons and alignment lengths are explained later 
in this document. My struggles with the expression data for these pseudogenes are described
throughout my 2022 lab notebook, especially Jan. 14 and Feb. 18.


I did this work in a conda environment that includes exonerate, macse, and some
other programs, but I think those two are all you need.

To activate this environment, do:
module load python/3.9-anaconda-2021.11
source activate /global/cfs/cdirs/plantbox/hybridum/software/exoxo






############# Identify remnants of 'lost genes' and get dN/dS relative to diploid ortholog #############

Let's start with the simplest and the original analysis: picking up 
unannotated pseudogenes in Bhyb26 by aligning a diploid peptide 
to the appropriate candidate region where we think the pseudogene is.
This is on cori in the dir 
/global/cscratch1/sd/vstartag/Bhyb26_lost_genesNEWNEW/lost_genes_dipsvsBhyb26 .

Start with the following data files in your working directory:
Bdistachyon_556_v3.2.cds_primaryTranscriptOnly.fa
Bdistachyon_556_v3.2.gene.gff3
Bdistachyon_556_v3.2.protein_primaryTranscriptOnly.fa
Bdistachyon__v__Bhyb26D.tsv
BhybridumBhyb26v2.1.gene.gff3
Brachypodium_hybridum_var_Bhyb26.mainGenome.fasta
Bstacei_316_v1.1.cds_primaryTranscriptOnly.fa
Bstacei_316_v1.1.gene.gff3
Bstacei_316_v1.1.protein_primaryTranscriptOnly.fa
Bstacei__v__Bhyb26S.tsv
gffWithOgs.txt
lostgenes.txt
^Some of these files are from phytozome, some are from GENESPACE, and some
are from my script bin_orthogroups.R, which processes GENESPACE output.

You'll also need all the scripts in your working directory. Then do:
python3 genespace2intervalsNEW.py
bedtools sort -i missing_gene_intervals.bed > missing_gene_intervals.sorted.bed
bedtools getfasta -fi Brachypodium_hybridum_var_Bhyb26.mainGenome.fasta -bed missing_gene_intervals.sorted.bed -nameOnly -fo missing_gene_intervals.fa
python3 parseRegionsFromFasta.py
python3 parsePep.py
python3 make_slurm_aliases.py #confusing name because I'm not actually using these for slurm, it's a holdover from an earlier v. of pipeline
./runExonerate_p2g.sh
for fname in $( ls candidatePseudogeneRegions/*.aln.p2g ); do python3 extract_seq2_exonerate.py $fname ; done
python3 parseCDS.py
./runMACSE.sh
rm candidatePseudogeneRegions/*.cds_macse_AA.fasta #remove extra stuff we don't need
for fname in $( ls candidatePseudogeneRegions/*.cds_macse_NT.fasta ); do python3 changeStopsToAmbigCharsNEW.py $fname ; done
python3 makeCtlFileNEW.py
./run_yn00.sh
python3 extractdNdSNEW.py


Here's an overview of what those commands did:
I identified the genomic neighborhood where the missing gene 'should be' with my 
script genespace2intervalsNEW.py. This produces the file missing_gene_intervals.bed,
which contains the candidate region coordinates and the names of the corresponding 
diploid genes. Then I got all of those genomic regions from the Bhyb26 genome as one big fasta
file with bedtools. Next, I split the big fasta file into many little fasta files--one per candidate region.
I put these in a directory called candidatePseudogeneRegions. Then for each corresponding
progenitor ortholog, I extracted the primary peptide sequence from the 
*protein_primaryTranscriptOnly.fa file from Phytozome.
So now in candidatePseudogeneRegions we have one Bhyb26 gDNA fasta file and one diploid 
peptide fasta file for each "missing gene".

Then I ran a quick script to rename my files for easy submission to slurm
as a job array, and I create a file listing the old names and new names 
called slurm_aliases.txt. This ended up not being necessary, since exonerate was so 
fast that I didn't need to submit a job for it, so it's been kind of grandfathered 
in to this pipeline. However, it is handy and tidy to have a number for each pseudogene/syn. orthogroup.
Though note each analysis directory has its own slurm aliases that don't correspond to each other.
Next, I ran exonerate to align the diploid peptide to the Bhyb26 genomic region. 
It only takes a minute or two. Next, I pull from the exonerate output the Bhyb26 sequence that 
could be aligned to the diploid peptide. I then align the diploid CDS to this Bhyb26 DNA sequence.
This second alignment step ensures that we have an in-frame DNA-DNA alignment we can
feed to PAML/yn00. I found that MACSE's codon-aware alignment was much better than exonerate's,
even though exonerate worked beautifully on the peptide to genome alignment.
It takes like half an hour to an hour to run MACSE on all the pseudogenes.
Next, I cleaned up the MACSE alignments a bit because any stop codon (TGA, TAA, or TAG)
will break yn00. Some of these stop codons are real premature termination codons (PTCs)
in Bhyb26, and some of them are spurious from questionable alignments. See my lab notebook
from 10 Mar. 2022 for more on this. Finally, I ran yn00 to get dN/dS values.

To get PTC (premature termination codon) information from the putative pseudogenes, do:
for fname in $( ls candidatePseudogeneRegions/*.aln.p2g ); do python3 PTCs_from_exonerate.py $fname ; done
Info on the length of the alignment versus the length of the initial diploid peptide/CDS
is already produced by extract_seq2_exonerate.py, which produces a file called length_differences.txt.

##########################################################################################








############# Get dN/dS for random samples of fully conserved Bhyb26 genes, relative to diploid ortholog #############

Cool beans. So now that we have dN/dS values for Bhyb26 "missing genes" relative to their diploid
ortholog, we want to get well-conserved Bhyb26 relative to their diploid ortholog as a control. 
Since 464 genes is not a very large sample size, we will do the whole procedure 1,000 times.
So we are randomly sampling 464 fully conserved Bhyb26 genes and getting stats on those,
then we are repeating the sampling and analysis 999 more times.

control_genes_dipvsBhyb26:

You will need the following data files in your working dir:

Bdistachyon_556_v3.2.cds_primaryTranscriptOnly.fa
Bdistachyon_556_v3.2.gene.gff3
Bdistachyon_556_v3.2.protein_primaryTranscriptOnly.fa
Bdistachyon__v__Bhyb26D.tsv
BhybridumBhyb26v2.1.gene.gff3
Brachypodium_hybridum_var_Bhyb26.mainGenome.fasta
Brachypodium_hybridum_var_Bhyb26.mainGenome.fasta.fai
Bstacei_316_v1.1.cds_primaryTranscriptOnly.fa
Bstacei_316_v1.1.gene.gff3
Bstacei_316_v1.1.protein_primaryTranscriptOnly.fa
Bstacei__v__Bhyb26S.tsv
conservedgenes.txt # from my script bin_orthogroups.R


Here's the batch script to run 1,000 trials:

#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH -t 1:30:00
#SBATCH -q genepool_shared
#SBATCH --array=1-10 #do 1000
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mail-user=vstartaglio@lbl.gov
#SBATCH --mail-type=end

module load python/3.9-anaconda-2021.11
source activate /global/cfs/cdirs/plantbox/hybridum/software/exoxo

python3 choose_conserved_ctrlsNEW.py sample${SLURM_ARRAY_TASK_ID} 224 240
bedtools sort -i sample${SLURM_ARRAY_TASK_ID}/conserved_gene_intervals.bed > sample${SLURM_ARRAY_TASK_ID}/conserved_gene_intervals.sorted.bed
bedtools getfasta -fi Brachypodium_hybridum_var_Bhyb26.mainGenome.fasta -bed sample${SLURM_ARRAY_TASK_ID}/conserved_gene_intervals.sorted.bed -nameOnly -fo sample${SLURM_ARRAY_TASK_ID}/conserved_gene_intervals.fa
python3 parseRegionsFromFasta.conserved.py sample${SLURM_ARRAY_TASK_ID}
python3 parsePep.conserved.py sample${SLURM_ARRAY_TASK_ID}
./runExonerate_p2g.conserved.sh sample${SLURM_ARRAY_TASK_ID}
for fname in $( ls sample${SLURM_ARRAY_TASK_ID}/candidatePseudogeneRegions/*.aln.p2g ); do python3 extract_seq2_exonerate.py $fname ; done
python3 parseCDS.conserved.py sample${SLURM_ARRAY_TASK_ID}
./runMACSE.conserved.sh sample${SLURM_ARRAY_TASK_ID}
rm sample${SLURM_ARRAY_TASK_ID}/candidatePseudogeneRegions/*.cds_macse_AA.fasta 
for fname in $( ls sample${SLURM_ARRAY_TASK_ID}/candidatePseudogeneRegions/*.cds_macse_NT.fasta ); do python3 changeStopsToAmbigCharsNEW.py $fname ; done
python3 makeCtlFile.conservedNEW.py sample${SLURM_ARRAY_TASK_ID}
cd sample${SLURM_ARRAY_TASK_ID}
cp ../run_yn00.conservedNEW.sh .
./run_yn00.conservedNEW.sh 
cd ..
python3 extractdNdS.conserved.py sample${SLURM_ARRAY_TASK_ID}



Finally, to get all the PTC (premature termination codon) info in one file, do:
sbatch PTCs_from_exonerate.sh #takes about 33 hours! Could probably speed it up by submitting
#a handful of jobs covering e.g. 250 files.
for i in $(seq 1 1000); do python3 gather_control_PTCs.py sample${i}/PTCs_from_exonerate.txt ; done
Also do this to get the average difference between the initial peptide length
and the final alignment length, for each trial:
for i in $(seq 1 1000); do python3 gather_control_lengths.py sample${i}/length_differences.txt ; done
Here's a little python code to get the final average of all trials, so an average of averages:

res1 = 0
res2 = 0
counter = 0 # should reach 1000...
with open('PTC_summary.txt', 'r') as inF:
	read = inF.readlines()
	for n in range(1, len(read)):
		L = read[n].strip('\n').split('\t')
		res1 += int(L[1])
		res2 += int(L[2])
		counter += 1

res1 #total number of PTCs in diploid gene from alignments
res2 #total number of PTCs in polyploid gene from alignments
round(res1/counter, 3) #mean number of PTCs in diploid, per trial
round(res2/counter, 3) #mean number of PTCs in polyploid, per trial

Here's the same thing but for alignment lengths:
res = 0
counter = 0 # should reach 1000...
with open('length_summary.txt', 'r') as inF:
	read = inF.readlines()
	for n in range(1, len(read)):
		L = read[n].strip('\n').split('\t')
		res += float(L[1])
		counter += 1

round(res/counter, 3) #mean difference between diploid and alignment, all trials

for i in $(seq 1 500); do mv sample${i}/dNdSbyGene.txt sample${i}/dSbyGene.txt ; done


Next, get diploid orthologs of the lost genes, in brachy and rice.
We will align these CDSs to each other and get dN/dS values.
This will go into our "pseudogenization time" analysis.
/global/cscratch1/sd/vstartag/Bhyb26_lost_genesNEWNEW/lost_genes_dipsvsrice :
prep_rice.py
sbatch runMACSErice.sh
rm candidatePseudogeneRegions/*_macse_AA.fasta 
for fname in $( ls candidatePseudogeneRegions/*_macse_NT.fasta ); do python3 changeStopsToAmbigCharsNEW.py $fname ; done
python3 makeCtlFileNEW.py
./run_yn00.sh
python3 extractdNdSNEW.py

##########################################################################################









############# Get dN/dS relative to progenitor ortholog: ABR113 lost genes #############

Next: ABR113 vs. progenitor orthologs.
We are basically repeating the pipeline with ABR113 instead of Bhyb26.
Working in /global/cscratch1/sd/vstartag/Bhyb26_lost_genesNEWNEW/lost_genes_dipsvsABR113, starting w these files + scripts:
Bdistachyon_556_v3.2.cds_primaryTranscriptOnly.fa
Bdistachyon_556_v3.2.gene.gff3
Bdistachyon_556_v3.2.protein_primaryTranscriptOnly.fa
Bdistachyon__v__ABR113D.tsv
Bhybridum_463_v1.0.fa
Bhybridum_463_v1.1.gene.gff3
Bstacei_316_v1.1.cds_primaryTranscriptOnly.fa
Bstacei_316_v1.1.gene.gff3
Bstacei_316_v1.1.protein_primaryTranscriptOnly.fa
Bstacei__v__ABR113S.tsv
gffWithOgs.txt
lostgenes.ABR113.txt #made by my script bin_orthogroups.R

python3 genespace2intervalsNEW.ABR113.py
bedtools sort -i missing_gene_intervals.bed > missing_gene_intervals.sorted.bed
bedtools getfasta -fi Bhybridum_463_v1.0.fa -bed missing_gene_intervals.sorted.bed -nameOnly -fo missing_gene_intervals.fa
python3 parseRegionsFromFasta.py
python3 parsePep.py
python3 make_slurm_aliases.py #a confusing name because I'm not actually using these aliases for slurm, holdover from an earlier v. of pipeline
./runExonerate_p2g.sh
for fname in $( ls candidatePseudogeneRegions/*.aln.p2g ); do python3 extract_seq2_exonerate.py $fname ; done
python3 parseCDS.py
sbatch runMACSE.sh
rm candidatePseudogeneRegions/*_macse_AA.fasta 
for fname in $( ls candidatePseudogeneRegions/*_macse_NT.fasta ); do python3 changeStopsToAmbigCharsNEW.py $fname ; done
python3 makeCtlFileNEW.py
./run_yn00.sh
python3 extractdNdSNEW.py


Output of genespace2intervalsNEW.ABR113.py, for posterity:
Started with 251 orthogroups containing a putative lost gene.
251 orthogroups remained after selecting those with a single gene in the appropriate diploid.
251 orthogroups remained after requiring that there be ten flanking (or nearly flanking) diploid genes with a single ortholog in the appropriate Bhyb26 subgenome.
237 orthogroups remained after requiring that at least 4 of the 5 B. hybridum orthologs on either side of the diploid gene be within 200kb of each other.
236 orthogroups remained after requiring that the final upstream and downstream B. hybridum 'anchor' genes demarcating my candidate region be within 200kb of each other.
##########################################################################################












############# Get dN/dS relative to rice: Bhyb26 lost genes #############

Next, align Bhyb26 'lost genes' to their rice orthologs and get dN/dS values.
This will go into our "pseudogenization time" analysis.
Need gffWithOgs.txt and Osativa_323_v7.0.cds_primaryTranscriptOnly.fa in your working dir.
/global/cscratch1/sd/vstartag/Bhyb26_lost_genesNEWNEW/lost_genes_Bhyb26vsrice :
mkdir candidatePseudogeneRegions
cp ../lost_genes_dipsvsBhyb26/candidatePseudogeneRegions/*.aln.seq2out candidatePseudogeneRegions
cp ../lost_genes_dipsvsBhyb26/slurm_aliases.txt .
python3 prep_rice.Bhyb26.py
sbatch runMACSE.sh
rm candidatePseudogeneRegions/*_macse_AA.fasta 
for fname in $( ls candidatePseudogeneRegions/*_macse_NT.fasta ); do python3 changeStopsToAmbigCharsNEW.py $fname ; done
python3 makeCtlFileNEW.py
./run_yn00.sh
python3 extractdNdSNEW.py

#########################################################################














############# Get dN/dS relative to rice: ABR113 lost genes #############

Finally, align ABR113 'lost genes' to their rice orthologs and get dN/dS values.
This will go into our "pseudogenization time" analysis. This is basically identical
to the protocol for lost_genes_Bhyb26vsrice, and you can even use the same scripts.
Need gffWithOgs.txt and Osativa_323_v7.0.cds_primaryTranscriptOnly.fa in your working dir.
/global/cscratch1/sd/vstartag/Bhyb26_lost_genesNEWNEW/lost_genes_ABR113vsrice :
mkdir candidatePseudogeneRegions
cp ../lost_genes_dipsvsABR113/candidatePseudogeneRegions/*.aln.seq2out candidatePseudogeneRegions
cp ../lost_genes_dipsvsABR113/slurm_aliases.txt .
python3 prep_rice.Bhyb26.py
sbatch runMACSE.sh
rm candidatePseudogeneRegions/*_macse_AA.fasta 
for fname in $( ls candidatePseudogeneRegions/*_macse_NT.fasta ); do python3 changeStopsToAmbigCharsNEW.py $fname ; done
python3 makeCtlFileNEW.py
./run_yn00.sh
python3 extractdNdSNEW.py

#########################################################################



############# Get dN/dS relative to progenitor ortholog: ALL fully conserved Bhyb26 genes #############
What happens if I just do all fully conserved genes, instead of samplings?
Actually, I will sample one Bhyb26 gene per syntenic orthogroup, since 15,000
genes should be plenty.

I'm working in /global/cscratch1/sd/vstartag/Bhyb26_lost_genesNEWNEW/all_conserved_genes.
I'm running this in a couple of steps so I can make sure things are finishing properly.

Here's step 1:

#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH -t 2:45:00
#SBATCH -q genepool_shared
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mail-user=vstartaglio@lbl.gov
#SBATCH --mail-type=end

python3 choose_conserved_all.py
bedtools sort -i conserved_gene_intervals.bed > conserved_gene_intervals.sorted.bed
bedtools getfasta -fi Brachypodium_hybridum_var_Bhyb26.mainGenome.fasta -bed conserved_gene_intervals.sorted.bed -nameOnly -fo conserved_gene_intervals.fa
python3 parseRegionsFromFasta.conserved.all.py
./runExonerate_p2g.sh 
for i in $(seq 1 1000); do 
python3 extract_seq2_exonerate.py candidatePseudogeneRegions/${i}.aln.p2g
done
python3 parseCDS.py

The next part will be slower so let's do it in parallel chunks. Here's my script step2.sh:
#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH -t 2:20:00 
#SBATCH -q genepool_shared
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mail-user=vstartaglio@lbl.gov
#SBATCH --mail-type=end

module load python/3.9-anaconda-2021.11
source activate /global/cfs/cdirs/plantbox/hybridum/software/exoxo

for ind in $(seq $start $end); do
if test -f "candidatePseudogeneRegions/${ind}.aln.seq2out"; then
	macse -prog alignSequences -seq candidatePseudogeneRegions/${ind}.cds.fa -seq_lr candidatePseudogeneRegions/${ind}.aln.seq2out > MACSEoutput.txt
fi
done

^Run this on the command line like so:
sbatch --export=ALL,start=1,end=1000 step2.sh
I ran these in chunks of 1,000 genes like so:
sbatch --export=ALL,start=1001,end=2000 step2.sh
sbatch --export=ALL,start=2001,end=3000 step2_1.sh
sbatch --export=ALL,start=3001,end=4000 step2_2.sh
sbatch --export=ALL,start=4001,end=5000 step2.sh
sbatch --export=ALL,start=5001,end=6000 step2_1.sh
sbatch --export=ALL,start=6001,end=7000 step2_2.sh
sbatch --export=ALL,start=7001,end=8000 step2.sh
sbatch --export=ALL,start=8001,end=9000 step2_1.sh
sbatch --export=ALL,start=9001,end=10000 step2_2.sh
sbatch --export=ALL,start=10001,end=11000 step2.sh
sbatch --export=ALL,start=11001,end=12000 step2_1.sh
sbatch --export=ALL,start=12001,end=13000 step2_2.sh
sbatch --export=ALL,start=13001,end=14000 step2_3.sh
sbatch --export=ALL,start=14001,end=15167 step2_4.sh

Next, finish it off. Let's call this step3:

#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH -t 2:00:00
#SBATCH -q genepool_shared
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mail-user=vstartaglio@lbl.gov
#SBATCH --mail-type=end

module load python/3.9-anaconda-2021.11
source activate /global/cfs/cdirs/plantbox/hybridum/software/exoxo

rm candidatePseudogeneRegions/*_macse_AA.fasta 
for fname in $( ls candidatePseudogeneRegions/*_macse_NT.fasta ); do python3 changeStopsToAmbigCharsNEW.py $fname ; done
python3 makeCtlFileNEW.py
./run_yn00.sh
python3 extractdNdSNEW.py

#########################################################################




