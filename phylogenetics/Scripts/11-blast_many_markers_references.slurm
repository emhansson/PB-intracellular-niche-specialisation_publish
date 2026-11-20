#!/bin/bash --login

#SBATCH -p multicore
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --job-name=blast_all_markers_references
#SBATCH -o blast_all_markers_references_%j.out
#SBATCH --error=blast_all_markers_references_%j.err

module purge

# Load BLAST
module load apps/binapps/blast/2.13.0

###############################################################################
# DIRECTORIES
###############################################################################

PROJECT=PROJECT_FOLDER

ASSEMBLIES=$PROJECT/references/reference_plastomes

MARKERS=$PROJECT/markers/NC_071969

OUTDIR=$PROJECT/blast_all_markers_reference

DBDIR=$PROJECT/references/reference_blastdb

assembly_type=reference

mkdir -p $DBDIR

mkdir -p $OUTDIR

###############################################################################
# SUMMARY FILE
###############################################################################

SUMMARY=$OUTDIR/blast_marker_summary.tsv

echo -e "sample_id\tmarker\tassembly_type\tn_hits\ttop_contig\ttop_identity\ttop_length\ttop_bitscore\ttop_evalue" \
> $SUMMARY

###############################################################################
# LOOP THROUGH ASSEMBLIES
###############################################################################

   for assembly in $ASSEMBLIES/*.fasta
do

    sample=$(basename "$assembly" .fasta)

    dbbase=$DBDIR/${sample}

    mkdir -p $OUTDIR/$sample


    ###########################################################################
    # MAKE BLAST DATABASE
    ###########################################################################

       makeblastdb \
        -in "$assembly" \
        -dbtype nucl \
        -out "$dbbase"

    ###########################################################################
    # BLAST MARKERS
    ###########################################################################

    for markerfasta in $MARKERS/*.fasta
    do

        marker=$(basename $markerfasta .fasta)

        echo "  BLASTing $marker"

        outfile=$OUTDIR/$sample/${marker}_hits.tsv

        blastn \
            -query $markerfasta \
            -db $dbbase \
            -out $outfile \
            -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
            -max_target_seqs 10

        #######################################################################
        # ADD TOP HIT TO SUMMARY
        #######################################################################

        if [[ -s "$outfile" ]]
        then

            top_hit=$(sort -k12,12nr $outfile | head -n 1)

            n_hits=$(wc -l < $outfile)

            contig=$(echo "$top_hit" | cut -f2)
            identity=$(echo "$top_hit" | cut -f3)
            length=$(echo "$top_hit" | cut -f4)
            evalue=$(echo "$top_hit" | cut -f11)
            bitscore=$(echo "$top_hit" | cut -f12)

            echo -e "${sample}\t${marker}\t${assembly_type}\t${n_hits}\t${contig}\t${identity}\t${length}\t${bitscore}\t${evalue}" \
            >> $SUMMARY

        else

            echo -e "${sample}\t${marker}\t${assembly_type}\t0\tNA\tNA\tNA\tNA\tNA" \
            >> $SUMMARY

        fi

    done

done

echo "All marker BLASTs complete."
echo "Summary table:"
echo $SUMMARY