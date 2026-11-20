#!/bin/bash --login

#SBATCH -p multicore 
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --job-name=getorg_array
#SBATCH --array=1-157%10
#SBATCH -o getorg_array_%j.out
#SBATCH --error=getorg_array_%j.err

# This script makes the organelle assemblies using default settings

module purge
source /opt/apps/apps/binapps/conda/25.3.0/etc/profile.d/conda.sh
conda activate getorganelle_env

PROJECT=PROJECT_FOLDER
MANIFEST=$PROJECT/qc/fastqc_temp/samples.tsv
OUTDIR=$PROJECT/getorganelle/all_runs

mkdir -p $OUTDIR

#Pull the correct line from the manifest
LINE=$(tail -n +2 $MANIFEST | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Extract columns
SAMPLE_ID=$(echo "$LINE" | cut -f1)
R1=$(echo "$LINE" | cut -f2)
R2=$(echo "$LINE" | cut -f3)

echo "Running sample: $SAMPLE_ID"
echo "R1: $R1"
echo "R2: $R2"

get_organelle_from_reads.py \
    -1 $R1 \
    -2 $R2 \
    -o $OUTDIR/${SAMPLE_ID}_getorganelle \
    -F embplant_pt \
    -t 8

