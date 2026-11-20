## QC
Run 01-fastqc.slurm, then run 02-multiqc.slurm

Just eyeballing the report, it's hard to tell which samples are the best, so moving onto chloroplast assembly so I can use both QC data and assembly quality to judge which sample to keep

## Chloroplast assembly
Using GetOrganelle as the 186b assembly only yielded a partially assembled chloroplast
03-make_manifest.py extracted all the symlinks from qc/fastqc_temp to create samples.tsv of where each sample is with shorter sample names
There's a total of 157 samples
04-all_getorganelle.slurm was used to apply getorganelle across this array

## Putting together a spreadsheet for comparing the samples
Used 05-get_assembly_stats.py, 06-add_metadata_to_assembly_stats.py and 07-make_final_qc_assembly_stats_table.py to make a big table of all assemblies to get a sense of how well it's worked (then took this to R: 08-QCInspection.R).

Should probably have just had one script for all of this and not bothered with the intermediates, or piped it or something but hey ho. 

Note that there was a typo in the pond_info.tsv where a pond is recorded as RX10, this should be RIX10 and is corrected later in the R script for downstream plotting.

About half are excellent/likely useable. All but one of the remaining has a stupid number of contigs and total size, but a largest contig that is nice and long and matches the expected plastid size -- looks like assembly of the actual plastid has worked quite well actually, but also ended up with a bunch of other stuff in there/graph explosion.

  assembly_quality             n proportion
  <chr>                    <int>      <dbl>
1 excellent                   39     0.419 
2 likely_usable                7     0.0753
3 poor                         1     0.0108
4 possible_marker_recovery    46     0.495 

This rules out trying to do something with the whole plastids, so moving onto marker recovery and then filtering what to include based on how good the recovery looks.

Will try later if needed: optimising the assemblies a bit for the really rough assemblies
-R 10 (or even 5, controls extension rounds, can reduce runaway graph recruitment)
-k 85, 105 (can reduce graph fragmentation for messy data)
--max-reads 5e6 (or even 2e6, since chloroplast coverage is high it can inflate errors and recruit contaminants and create the mess we seemingly have here, so this is worth testing first)

(Didn't end up doing the above, the marker recovery was overall excellent regardless of assembly quality)

## Testing marker gene extraction
Using the M. simp ref (NC_071969, https://www.ncbi.nlm.nih.gov/nuccore/NC_071969.1?report=gbwithparts&log$=seqview) Guy used to guide the 186b reference assembly.

inspect_gene_names.py for finding all the named genes in the annotation (just in case of any weird naming, but also useful to see what else there is)
extract_markers.py to actually extract the markers I want from the reference

Used 09-blast_many_markers.slurm to blast markers for all samples after first testing a subset picked based on assembly quality.

This was really successful, and shows that even the really bad assemblies actually seem to recover markers really well -- where there are duplicates, they either agree or one is of markedly worse quality overall.

Inspection of the marker recovery (based on marker_recovery_summary.tsv output) is in 10-MarkerRecovery.R

From this, a final list of best samples and a final list of core markers was created.

11-blast_many_marker_references.slurm was used to blast reference assemblies after the creation of the first exploratory tree (tree1) that only had field strains in it.
Based on recovery success, tree2 was constructed with Micractinium, Chlorella and Parachlorella references included.
Tree3 also included Chlamydomonas, Dunliella and Choricystis as outgroups as well as our 186b reference.

# Marker extraction
Actual marker extraction done with 12-extract_marker_fastas.slurm for field strains.
The markers for the references were appended to the extracted marker fastas with 12-add_more_extracted_marker_fastas.slurm

# Aligning, trimming, concatenation of markers
Aligning done with 14-align_markers.slurm
Trimming done with 15-trim_markers.slurm (this step was skipped for tree1 as spot checks of alignments looked very healthy, included for trees2 and 3, but not tree4. Trimming doesn't seem to have had any effect on topology)

Tom advised not trimming, so tree4 (final tree used in publication) is on the untrimmed alignments

For concatenation:

Activate amas_env
$ conda activate amas_env

# Run the below code
$ python -m amas.AMAS concat -i trimmed_alignments/*_aligned_trimmed.fasta -f fasta -d dna -p concatenation/partitions.txt -t concatenation/concatenated_alignment.fasta

Need to fix the formatting of the partitions file for iqtree to work though:
awk '{print "DNA, "$0}' partitions.txt > iqtree_partitions.txt

## Make phylogeny
Run 16-iqtree.slurm

## Marker trees
17-marker_trees.slurm makes a tree for each marker, this is used in 17-Trees.R for calculating evolutionary rates

## Plotting and calculation of marker rates
17-Trees.R for handling of trees and analysis of marker variation



