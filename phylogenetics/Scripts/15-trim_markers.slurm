#!/bin/bash --login

#SBATCH -p multicore
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --job-name=trimal

module purge
source /opt/apps/apps/binapps/conda/25.3.0/etc/profile.d/conda.sh
conda activate tree_env

PROJECT=PROJECT_FOLDER

ALIGNDIR=$PROJECT/aligned_markers
OUTDIR=$PROJECT/trimmed_alignments

mkdir -p $OUTDIR

for aln in $ALIGNDIR/*_aligned.fasta
do

    base=$(basename "$aln" .fasta)

    trimal \
        -in "$aln" \
        -out "$OUTDIR/${base}_trimmed.fasta" \
        -automated1

done