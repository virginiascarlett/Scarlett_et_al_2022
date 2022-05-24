#!/bin/bash
FILENUM=`wc -l ${1}/slurm_aliases.txt`
FN=`echo $FILENUM | sed "s/[^0-9]//g"`
for ind in $(seq 1 $FN); do
if test -f "${1}/candidatePseudogeneRegions/${ind}.aln.seq2out"; then
	macse -prog alignSequences -seq ${1}/candidatePseudogeneRegions/${ind}.cds.fa -seq_lr ${1}/candidatePseudogeneRegions/${ind}.aln.seq2out #> MACSEoutput.txt
fi
done

#STARTED WRITING AS A BATCH SCRIPT BUT ACTUALLY I THINK IT'S FASTER TO JUST RUN INTERACTIVELY