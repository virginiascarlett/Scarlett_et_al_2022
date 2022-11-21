# Scarlett_et_al_2022
Here you will find the scripts that were used in Scarlett et al. "Multiple origins, one evolutionary trajectory: gradual evolution characterizes distinct lineages of allotetraploid Brachypodium" Genetics, 2022, iyac146, https://doi.org/10.1093/genetics/iyac146. This repository includes the scripts that produced the figures, and where possible, also includes the analysis pipelines. See the data availability section of the paper for the links to the raw data themselves.

Contents:

fig1: Code to run GENESPACE

fig2: Code to produce GENESPACE dotplots

fig3: Code to identify putative pseudogenes and analyze their characteristics. There are two directories: one with the scripts I actually used, and one with tidier scripts and better documentation. If you wish to adapt this protocol for your own needs, use the tidied version. This protocol was designed for tetraploids with well-defined subgenomes and either one or two progenitors.

fig4: Code to map RNA-seq reads with BBmap, analyze the counts, calculate TPMs, and plot bar charts. See the comments at the top of countsToTPMbasicNEW.py for an example HTSeq call as well. 

fig5: A very messy directory containing code for (1) getting binned genome stats (e.g. gene density) and plotting them with circos, and (2) running my TE annotation pipeline and plotting some of the results.
(1) The circos/ subfolder is not too messy. getGeneDensity.R should be pretty reusable. The other get...R scripts are hastily modified versions of getGeneDensity.R.
I never tidied this for reproducibility since I think EDTA (https://github.com/oushujun/EDTA) is similar to my TE annotation pipeline, but better. I do not recommend you try to run my pipeline, since it wasn't designed as a stable reusable tool, and it was designed around our particular HPC cluster. Nevertheless, all the scripts I used are here. The README is more of a notebook than a tutorial, documenting what I did, but not how you can reproduce it. The scripts directly under fig5/ are for plotting the bar charts and the insertion time scatterplot. 

