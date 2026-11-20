#!/bin/bash --login

#SBATCH -p multicore
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --job-name=marker_trees
#SBATCH --array=1-19
#SBATCH -o marker_trees_+%j.out
#SBATCH --error=marker_trees_%j.err

module purge

source /opt/apps/apps/binapps/conda/25.3.0/etc/profile.d/conda.sh
conda activate tree_env

which iqtree

PROJECT=PROJECT_FOLDER

OUTDIR=$PROJECT/marker_trees

mkdir -p $OUTDIR

relpath=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$PROJECT/marker_list.txt")
aln="$PROJECT/$relpath"

echo "Task: $SLURM_ARRAY_TASK_ID"
echo "Alignment: $aln"

ls -l "$aln"

marker=$(basename "$aln" _aligned.fasta)

iqtree \
    -s "$aln" \
    -m GTR \
    -nt "$SLURM_CPUS_PER_TASK" \
    -pre $OUTDIR/"$marker"