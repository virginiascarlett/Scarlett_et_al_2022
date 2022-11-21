#!/bin/bash
#SBATCH -N 1
#SBATCH --cpus-per-task=64
#SBATCH --ntasks=1
#SBATCH --mem=118G
#SBATCH -t 72:00:00
#SBATCH -q genepool
#SBATCH -A plant
#SBATCH -C haswell
#SBATCH --mail-user=vstartaglio@lbl.gov
#SBATCH --mail-type=end
#SBATCH --array=1-26
#SBATCH --error=mcclint.%A.%a.err
#SBATCH --output=mcclint.%A.%a.out

CONDA_BASE=$(conda info --base)
source ${CONDA_BASE}/etc/profile.d/conda.sh
conda activate /global/homes/j/jespi/.conda/envs/mcclintock

INPUTDIR=/global/projectb/scratch/vstartag/jacob_internship/McClintock_inputs
WORKDIR=/global/projectb/scratch/vstartag/jacob_internship/Bhybridum
FASTQDIR=/global/projectb/scratch/vstartag/jacob_internship/jacob_fastqs

python3 /global/cscratch1/sd/jespi/mcclintock/mcclintock.py \
	-r $INPUTDIR/Bhybridum_463_v1.0.unmasked.fa \
	-c $INPUTDIR/TElib.ABR113.clean.fa \
	-1 $FASTQDIR/mcclintock${SLURM_ARRAY_TASK_ID}_R1.fastq.gz \
	-2 $FASTQDIR/mcclintock${SLURM_ARRAY_TASK_ID}_R2.fastq.gz \
	-g $INPUTDIR/ABR113_RepeatMasker_TElib.clean.gff \
	-t $INPUTDIR/ABR113_TElib_McClintockfams.tsv \
	-m temp2 \
	-p 64 \
	-o $WORKDIR/results_${SLURM_ARRAY_TASK_ID}
