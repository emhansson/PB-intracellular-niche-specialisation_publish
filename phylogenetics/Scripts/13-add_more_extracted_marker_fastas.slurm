#!/bin/bash --login

#SBATCH -p multicore
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --job-name=add_more_extract_markers
#SBATCH -o add_more_extract_markers_%j.out
#SBATCH --error=add_more_extract_markers_%j.err

# THIS WILL ADD TO EXISTING FILES SO BE CAREFUL TO ONLY RUN ON ONES THAT HAVEN'T BEEN RUN BEFORE
# Just doing this to not have to recompute all strains when I just want to add the reference strains...

module purge

# Load tools
module load apps/gcc/samtools/1.21

###############################################################################
# DIRECTORIES
###############################################################################

PROJECT=PROJECT_FOLDER

REFDIR=$PROJECT/references/reference_plastomes
BLASTDIR=$PROJECT/blast_all_markers_reference
OUTDIR=$PROJECT/extracted_markers

mkdir -p $OUTDIR

###############################################################################
# INPUT FILE
###############################################################################

MARKER_LIST=$PROJECT/marker_lists/core_markers.txt

###############################################################################
# LOOP THROUGH REFERENCE FASTAS
###############################################################################

for assembly in \
    $REFDIR/Micractinium*.fasta \
    $REFDIR/Chlorella*.fasta \
    $REFDIR/Parachlorella*.fasta \
    $REFDIR/Chlamydomonas*.fasta \
    $REFDIR/Choricystis*.fasta \
    $REFDIR/Dunaliella*.fasta \
    $REFDIR/m_conductrix*.fasta 
do

    sample=$(basename "$assembly" .fasta)

    echo "Processing $sample"

    ###########################################################################
    # INDEX FASTA
    ###########################################################################

    samtools faidx "$assembly"

    ###########################################################################
    # LOOP THROUGH MARKERS
    ###########################################################################

    while read marker
    do

        blastfile="$BLASTDIR/${sample}/${marker}_hits.tsv"

        echo "Checking blast file:"
echo "$blastfile"

ls -lh "$blastfile"

        if [[ ! -s "$blastfile" ]]
        then
            echo "No hits for $sample $marker"
            continue
        fi

        #######################################################################
        # GET TOP HIT
        #######################################################################

        top_hit=$(sort -k12,12nr "$blastfile" | head -n 1)

        contig=$(echo "$top_hit" | cut -f2)
        start=$(echo "$top_hit" | cut -f9)
        end=$(echo "$top_hit" | cut -f10)

        outfile="$OUTDIR/${marker}.fasta"

        #######################################################################
        # EXTRACT SEQUENCE
        #######################################################################

        if (( start <= end ))
        then

            samtools faidx \
                "$assembly" \
                "${contig}:${start}-${end}" \
            | sed "1s/>.*/>${sample}/" \
            >> "$outfile"

        else

            samtools faidx \
                -i \
                "$assembly" \
                "${contig}:${end}-${start}" \
            | sed "1s/>.*/>${sample}/" \
            >> "$outfile"

        fi

        echo "Extracted $marker from $sample"

    done < "$MARKER_LIST"

done

echo "All reference marker extraction complete."