To create my circos plots, I used the following files:
-Karyotype files made manually (see circos documentation)
-TE density files, e.g. TEdensity.Bhyb26.200kb.tsv, which look like this:
BhD1    1       200000  0.065
BhD1    200001  400000  0.184
BhD1    400001  600000  0.115
BhD1    600001  800000  0.155
BhD1    800001  1000000 0.135
These come from my script getTEdensity.R, which requires:
  allfragments.classified_ABR113 #produced from step6 of my TE pipeline
  allfragments.classified_Bhyb26 #produced from step6 of my TE pipeline
-TE diversity files, e.g. TEdiversity.Bhyb26.200kb.circos.tsv, which look like this:
BhD1    1       200000  40
BhD1    200001  400000  53
BhD1    400001  600000  51
BhD1    600001  800000  44
BhD1    800001  1000000 46
BhD1    1000001 1200000 47
BhD1    1200001 1400000 37
BhD1    1400001 1600000 49
These come from my script getTEdiversity.R, which requires:
  allfragments.classified_ABR113 
  allfragments.classified_Bhyb26 
  ABR113.chrom.sizes #A two-column file with the chromosome in one column and its size, in bp, in the other
  Bhyb26.chrom.sizes #A two-column file with the chromosome in one column and its size, in bp, in the other
Note that TE diversity is measured as number of families per bin (in this case, 200kb). 
-Gene density files, which come from my script getGeneDensity.R
-Centromere and pericentromere files (see below for more info)
-TE density files for subgenome-specific TEs only. 
These were produced with get_subg_specific_TE_density.R, which requires:
    allfragments.classified_ABR113 
    allfragments.classified_Bhyb26 
    ABR113.chrom.sizes 
    Bhyb26.chrom.sizes 
    subg_specificity_by_family.txt_ABR113 #produced from step6 of my TE pipeline
    subg_specificity_by_family.txt_Bhyb26 #produced from step6 of my TE pipeline
Note that these files illustrate the distribution of TE copies that come from TE families with
more than 90% of TE copies on one subgenome, and 5 or more TE copies total.



To get a track under the ideogram illustrating the centromeres and pericentromeres, I made two overlapping histogram tracks 
with 0s and 1s as values. I used three scripts to prepare the input files. The first two can be run in either order, but 
make_01_histo_circos.R must be run last. chromosomes_to_bins.R will parse the karyotype file into a file with bins, e.g.:
BhD1,1,100000
BhD1,100001,200000
BhD1,200001,300000
BhD1,300001,400000
Note that currently, the bin size is hard-coded into this script. 

prepare_centro_files_circos.R. will make non-overlapping binned centromere and pericentromere files. 
As input this script requires e.g. ABR113centromeres.txt and Bhyb26centromeres.txt. These files were
produced manually from visual assessment of gene and TE density. They delineate the centromere and 
pericentromere locations in Mb, and look like this:
BhD1    27      49      pericentromere
BhD1    36      41      centromere
BhD2    22      36.5    pericentromere
BhD2    29      32      centromere
BhD3    18      36      pericentromere
BhD3    24      29.5    centromere
The output of prepare_centro_files_circos.R, e.g. Bhyb26centromeres.circos.txt, looks like this:
BhD1 36000000 41000000
BhD2 29000000 32000000
BhD3 24000000 29500000

make_01histo_circos.R will produce the files that I actually give to circos. Each one looks like:
BhD1 1 100000 0
BhD1 100001 200000 0
BhD1 200001 300000 0
BhD1 300001 400000 0
BhD1 400001 500000 0
etc., but Bhyb26centromeres.circos.final.txt will have 1's in the centromere bins and Bhyb26pericentromeres.circos.final.txt
will have 1's in the pericentromere bins.

I installed circos using a conda environment on the cori cluster.
I had some trouble installing circos, and I had to follow this post:
https://groups.google.com/g/circos-data-visualization/c/Guh7fjX8xHI
The way I created the environment was:
conda create --prefix /global/cfs/cdirs/plantbox/hybridum/software/bigtop circos
conda install -c conda-forge libgd=2.2.4
conda install -c bioconda perl-gd=2.56
