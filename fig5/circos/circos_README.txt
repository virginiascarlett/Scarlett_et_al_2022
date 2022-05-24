To create my circos plots, you need the following files:
allfragments.classified_ABR113 #produced from step6 of my TE pipeline
allfragments.classified_Bhyb26 #produced from step6 of my TE pipeline
ABR113.chrom.sizes #produced as part of my Hi-C pipeline
Bhyb26.chrom.sizes #produced as part of my Hi-C pipeline
subg_specificity_by_family.txt_ABR113 #produced from step6 of my TE pipeline
subg_specificity_by_family.txt_Bhyb26 #produced from step6 of my TE pipeline
ABR113centromeres.txt #produced manually and painstakingly from a combination of TE density and Hi-C
Bhyb26centromeres.txt #produced manually and painstakingly from TE density files

TE density files come from my script getTEdensity.R.
(On my personal laptop, that's in TEs/TE_pipeline/step6/getTEdensity.R)
TE diversity files come from my script getTEdiversity.R. TE diversity is measured as
number of families per bin (in this case, 200kb). (Find this script in densities/)
Subgenome-specific TE density files came from my script get_subg_specific_TE_density.R. 
These files illustrate the distribution of TE copies that come from TE families with
more than 90% of TE copies on one subgenome, and 5 or more TE copies total.

On my personal laptop (/github) there is a script in synteny/ called prepare_centro_files_circos.R.
That will convert the *.chrom.sizes files to a format circos can use.

To load environment, do:
module load python/3.9-anaconda-2021.11
source activate /global/cfs/cdirs/plantbox/hybridum/software/bigtop
To run circos, do:
/global/cfs/cdirs/plantbox/hybridum/software/bigtop/bin/circos

For posterity, I had some trouble installing circos, and I had to follow this post:
https://groups.google.com/g/circos-data-visualization/c/Guh7fjX8xHI
So the way I created the environment was:
conda create --prefix /global/cfs/cdirs/plantbox/hybridum/software/bigtop circos
conda install -c conda-forge libgd=2.2.4
conda install -c bioconda perl-gd=2.56

The final version that I ended up using was NOT the master annotations but the original,
genome-specific annotations (see Mar. 10, 2022 in my notebook).
