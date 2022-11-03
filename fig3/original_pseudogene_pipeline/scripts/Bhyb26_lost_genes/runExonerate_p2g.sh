#!/bin/bash
FILENUM=`wc -l slurm_aliases.txt`
FN=`echo $FILENUM | sed "s/[^0-9]//g"`
cd candidatePseudogeneRegions/
for ind in $(seq 1 $FN)
	do 
		exonerate --model protein2genome -n 1 --subopt False --showvulgar FALSE --exhaustive n ${ind}.pep.fa ${ind}.region.fa --ryo ">%qi %ql %qal\n%qas\n>%ti %tcl\n%tcs\n" > ${ind}.aln.p2g
	done
