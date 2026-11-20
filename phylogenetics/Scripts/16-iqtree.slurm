#!/bin/bash --login

#SBATCH -p multicore
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --job-name=plastid_tree
#SBATCH -o plastid_tree_%j.out
#SBATCH --error=plastid_tree_%j.err

module purge

source /opt/apps/apps/binapps/conda/25.3.0/etc/profile.d/conda.sh
conda activate tree_env

which iqtree

PROJECT=PROJECT_FOLDER

OUTDIR=$PROJECT/tree4_untrimmed

mkdir -p $OUTDIR

iqtree \
    -s $PROJECT/untrimmed_concatenation/concatenated_alignment.fasta \
    -p $PROJECT/untrimmed_concatenation/iqtree_partitions.txt \
    -m MFP+MERGE \
    -B 1000 \
    -T AUTO \
    -pre $OUTDIR/plastid_tree