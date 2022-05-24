#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=16
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH -t 5:00:00
#SBATCH -q genepool_shared
#SBATCH --array=3
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mail-user=vstartaglio@lbl.gov
#SBATCH --mail-type=end

module load python/3.9-anaconda-2021.11
source activate /global/cfs/cdirs/plantbox/hybridum/software/makeitcount
#htseq-count -r pos -f bam -s yes -t gene -i ID ${SLURM_ARRAY_TASK_ID}.bbmap_sorted.bam BhybridumBhyb26v2.1.gene_exons.gff3 > ${SLURM_ARRAY_TASK_ID}_HTScounts.txt
htseq-count -r pos -f bam -s yes -t gene -i ID ${SLURM_ARRAY_TASK_ID}.bbmap_sorted.bam Bhyb26gff_w_pseudogenes.gff3 > ${SLURM_ARRAY_TASK_ID}_HTScounts.pseudos.txt
