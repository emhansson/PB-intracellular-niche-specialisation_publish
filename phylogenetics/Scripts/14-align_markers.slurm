#!/bin/bash --login

#SBATCH -p multicore
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --job-name=mafft_markers
#SBATCH -o mafft_markers_%j.out
#SBATCH --error=mafft_markers_%j.err

module purge

# Load MAFFT
module load apps/intel-19.1/mafft/7.505

###############################################################################
# DIRECTORIES
###############################################################################

PROJECT=PROJECT_FOLDER

INDIR=$PROJECT/extracted_markers
OUTDIR=$PROJECT/aligned_markers

mkdir -p $OUTDIR

###############################################################################
# ALIGN ALL MARKERS
###############################################################################

for fasta in $INDIR/*.fasta
do

    marker=$(basename $fasta .fasta)

    echo "Aligning $marker"

    mafft \
        --auto \
        --thread 4 \
        "$fasta" \
        > "$OUTDIR/${marker}_aligned.fasta"

done

echo "All alignments complete."