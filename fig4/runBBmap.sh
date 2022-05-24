#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --ntasks=1
#SBATCH --mem=118G
#SBATCH -t 5:00:00
#SBATCH -q genepool
#SBATCH --array=1
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mail-user=vstartaglio@lbl.gov
#SBATCH --mail-type=end

module load python/3.7-anaconda-2019.10
source activate /global/projectb/sandbox/plant/hybridum/software/staycoolTEMP
STATS=mapping_stats

/global/cfs/cdirs/plantbox/hybridum/software/bbmap/bbmap.sh \
	in=${SLURM_ARRAY_TASK_ID}.fastq.gz \
	ref=Brachypodium_hybridum_var_Bhyb26.mainGenome.fasta \
	out=${SLURM_ARRAY_TASK_ID}.bbmap.sam \
	maxindel=100000 \
	bs=${SLURM_ARRAY_TASK_ID}.sam2bam.sh \
	refstats=$STATS/${SLURM_ARRAY_TASK_ID}.bbmap_refstats.txt \
	minid=0.9 \
	ambig=toss \
	covstats=$STATS/${SLURM_ARRAY_TASK_ID}.covstats.info \
	covhist=$STATS/${SLURM_ARRAY_TASK_ID}.covhist.plt \
	basecov=$STATS/${SLURM_ARRAY_TASK_ID}.basecov.info #; sh ${SLURM_ARRAY_TASK_ID}.sam2bam.sh
