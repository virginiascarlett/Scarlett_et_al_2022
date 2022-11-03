#!/bin/bash
filenum=`wc -l $1/slurm_aliases.txt`
stringarray=($filenum)
fn=${stringarray[0]}
#FN=`echo $FILENUM | sed "s/[^0-9]//g"`
cd $1/candidatePseudogeneRegions/
for ind in $(seq 1 $fn)
	do 
		exonerate --model protein2genome -n 1 --subopt False --showvulgar FALSE --exhaustive n ${ind}.pep.fa ${ind}.region.fa --ryo ">%qi %ql %qal\n%qas\n>%ti %tcl\n%tcs\n" > ${ind}.aln.p2g
	done
