# Scarlett_et_al_2022
Here you will find the scripts that were used in: 

Scarlett et al. "Multiple origins, one evolutionary trajectory: gradual evolution characterizes distinct lineages of allotetraploid Brachypodium" Genetics, 2022, iyac146, https://doi.org/10.1093/genetics/iyac146. 

This repository includes the scripts that produced the figures, and where possible, also includes the analysis pipelines. See the data availability section of the paper for the links to the raw data themselves.

Contents:

fig1: Code to run GENESPACE

fig2: Code to produce GENESPACE dotplots

fig3: Code to identify putative pseudogenes and analyze their characteristics. There are two directories: one with the scripts I actually used, and one with tidier scripts and better documentation. If you wish to adapt this protocol for your own needs, use the tidied version. This protocol was designed for tetraploids with well-defined subgenomes and either one or two progenitors.

fig4: Code to map RNA-seq reads with BBmap, analyze the counts, calculate TPMs, and plot bar charts. See table1/ for an example HTSeq call as well. 

fig5: A pretty messy directory containing code for (1) getting binned genome stats (e.g. gene density) and plotting them with circos, and (2) running my TE annotation pipeline and plotting some of the results. The scripts directly under fig5/ are for plotting the bar charts and the insertion time scatterplot. 

  (1) The circos/ subfolder has my scripts for producing circos input files. Confusingly, the script to produce the bar plot for fig 5 panel c is in here, at the bottom of get_subg_specific_TE_density.R. getGeneDensity.R is fairly tidy and reusable. Unfortunately, I was locked out of my NERSC/LBNL account before I could save a copy of my circos.conf files, so I don't have those to share.
  
  (2) The TE_annotation/ subfolder is messy. I never tidied this for reproducibility since EDTA (https://github.com/oushujun/EDTA) is similar to my TE annotation pipeline, but better. I do not recommend you try to run my pipeline, since it wasn't designed as a stable reusable tool, and it was designed around our particular HPC cluster. Nevertheless, all the scripts I used are here. The README is more of a notebook than a tutorial. It documents what I did, but not how you can reproduce it. 

suppfigS2: Code for visualizing the locations of putative pseudogenes on a Bhyb26 karyogram, and a histogram of the TPM values for conserved genes and putative pseudogenes. Note that the code that produced Sup. Fig. S2, panel c is at the bottom of fig5/hybridum_genes_distance_to_TEs.R. Note also that to calculate the TPMs for the putative pseudogenes, I just added my putative pseudogenes to the gff from Phytozome, and then I did a standard RNA-seq workflow with this modified gff. I don't think my script is reproducible, since the script that produced the modified gff, table1/intervals2gff.py, relies on an old file called pesudogeneIntervals.txt that my new pseudogene pipeline no longer produces. To reproduce this, the user will probably need to modify fig3/pseudogene_pipeline_tidied_for_reuse/scripts/genespace2intervals.py to get the actual locations in the genome of all the exons of the putative pseudogenes. 

suppfigS3: Code for running TEMP2 using the McClintock wrapper

table1: Code for getting TPMs for pseudogenes (see note under suppfigS2), and some general RNA-seq scripts. This is table 1 in the text, not supplementary table S1. Several of the columns in that supplementary table were calculated using fig5/hybridum_genes_distance_to_TEs.R.
