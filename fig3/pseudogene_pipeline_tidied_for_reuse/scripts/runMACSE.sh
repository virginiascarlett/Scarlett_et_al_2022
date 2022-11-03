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

FILENUM=`wc -l aliases.txt`
FN=`echo $FILENUM | sed "s/[^0-9]//g"`
for ind in $(seq 1 $FN); do
if test -f "candidatePseudogeneRegions/${ind}.aln.seq2out"; then
	macse -prog alignSequences -seq candidatePseudogeneRegions/${ind}.cds.fa -seq_lr candidatePseudogeneRegions/${ind}.aln.seq2out > MACSEoutput.txt
fi
done

#This is faster than submitting hundreds of 2-minute jobs