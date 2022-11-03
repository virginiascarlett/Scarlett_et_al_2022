#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH -t 33:00:00
#SBATCH -q genepool_shared
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mail-user=vstartaglio@lbl.gov
#SBATCH --mail-type=end

for i in $(seq 11 1000); do
	for fname in $( ls sample${i}/candidatePseudogeneRegions/*.aln.p2g ); do python3 PTCs_from_exonerate.py $fname ; done
done