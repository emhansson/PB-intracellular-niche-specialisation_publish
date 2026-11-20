#!/bin/bash --login

#SBATCH -p multicore
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --job-name=extract_markers
#SBATCH -o extract_markers_%j.out
#SBATCH --error=extract_markers_%j.err

module purge

# Load tools
module load apps/gcc/samtools/1.21

###############################################################################
# DIRECTORIES
###############################################################################

PROJECT=PROJECT_FOLDER

ASSEMBLIES=$PROJECT/getorganelle/all_runs
BLASTDIR=$PROJECT/blast_all_markers
OUTDIR=$PROJECT/extracted_markers

mkdir -p $OUTDIR

###############################################################################
# INPUT FILES
###############################################################################

MARKER_LIST=$PROJECT/marker_lists/core_markers.txt
SAMPLE_LIST=$PROJECT/sample_lists/best_samples.txt

###############################################################################
# LOOP THROUGH SAMPLES
###############################################################################

while read sample
do

    echo "Processing $sample"

    ###########################################################################
    # FIND ASSEMBLY FASTA
    ###########################################################################

    if [[ -f "$ASSEMBLIES/${sample}/embplant_pt.K115.complete.graph1.1.path_sequence.fasta" ]]
    then
        assembly="$ASSEMBLIES/${sample}/embplant_pt.K115.complete.graph1.1.path_sequence.fasta"

    elif [[ -f "$ASSEMBLIES/${sample}/embplant_pt.K115.scaffolds.graph1.1.path_sequence.fasta" ]]
    then
        assembly="$ASSEMBLIES/${sample}/embplant_pt.K115.scaffolds.graph1.1.path_sequence.fasta"

    else
        echo "No assembly found for $sample"
        continue
    fi

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

done < "$SAMPLE_LIST"

echo "All marker extraction complete."