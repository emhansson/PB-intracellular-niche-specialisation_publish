This FieldStrains_README.txt file was generated on 2026-11-12 by Erika M Hansson

GENERAL INFORMATION

1. Title of Dataset: Divergence between intracellular and extracellular niches drives endosymbiont evolutionary entrapment in Paramecium bursaria

2. Author Information
	A. Principal Investigator Contact Information
		Name: Michael A Brockhurst
		Institution: The University of Manchester
		Address: Division of Evolution, Infection and Genomics, 
		University of Manchester, Manchester M13 9PT, UK
		Email: michael.brockhurst@manchester.ac.uk

	B. Lead Author 1 Contact Information
		Name: Irma Vitonytė
		Institution: The University of Manchester
		Address: Division of Evolution, Infection and Genomics, 
		University of Manchester, Manchester M13 9PT, UK
		Email: irma.vitonyte@manchester.ac.uk

	C. Lead Author 2 and Author Responsible for the Code Contact Information
		Name: Erika M Hansson
		Institution: The University of Manchester
		Address: Division of Evolution, Infection and Genomics, 
		University of Manchester, Manchester M13 9PT, UK
		Email: erika.hansson@manchester.ac.uk

3. Date of data collection: Strains were collected from the field in summers 2021 and 2022. Characterisation of strains and DNA extraction occurred 2022--2025.

4. Geographic location of data collection: Greater Manchester area, United Kingdom. Exact locations of ponds available in the associated manuscript supplemental information.

5. Information about funding sources that supported the collection of the data: NERC (NE/V000128/1), BBSRC (BB/X016439/1).


SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: CC-BY-4.0

2. Recommended citation for this dataset: 

Irma Vitonytė, Erika M. Hansson, Daniel Malumphy-Montesdeoca, Guy Leonard, Fiona R. Savory, Luke T. Dunning, Thomas A. Richards, Andrew P. Beckerman, Duncan D. Cameron , Michael A. Brockhurst. (2026). Intracellular niche specialisation drives evolutionary entrapment of endosymbiotic algae in Paramecium bursaria. 

DATA & FILE OVERVIEW

1. File List: 

\phenotypes -- Everything related to analysis of phenotypic data

FieldIsolates_publish.Rproj -- RProject file for data analysis

	\Scripts
		\01-EnvData.R -- Summarising of pond chemistry data
		\02-NData.R -- Analysis of growth on different nitrogen sources
		\03-SugarData.R -- Analysis of sugar release in different pHs
		\04-pHGRData.R -- Analysis of effect of pH on growth, and whether this 					correlates with sugar release
		\05-FluorimetryData.R -- Analysis of Fv/Fm
		\06-MultivariateModelling.R -- Putting all the data together for a multivariate analysis
		\07-Plots.R -- Code for all the plots in the manuscript collected in one place

	\Data
		\df_tidy_env.csv -- Information about strain origins
			pond_id: originating pond
			pond_area: general area of pond
			strain: strain identifier
			lifestyle: free-living (FL) vs. symbiotic (S)
			species_id: strain species
			no3: nitrate levels measured in situ
			tds: tds levels measured in situ
			conductivity: conductivity measured in situ
			ph: pH measured measured in situ

		\df_tidy_n.csv
			strain: strain identifier
			lifestyle: free-living (FL) vs. symbiotic (S)
			tech_rep: technical replicate
			n_source: nitrogen source present in growth medium
			n_type: inorganic vs. organic
			reading: cell density reading (OD750)
			day: day of measurement (0 or 7)

		\df_tidy_sugar.csv
			strain: strain identifier
			lifestyle: free-living (FL) vs. symbiotic (S)
			tech_rep: technical replicate
			ph: 5.5 vs. 7.6
			sugar: Maltose, glucose or trehalose
			release: sugar release amount in ppm

		\df_tidy_pH.csv
			strain: strain identifier
			lifestyle: free-living (FL) vs. symbiotic (S)
			tech_rep: technical replicate
			ph: 5.5 vs. 7.6
			od750: cell density reading (OD750)
			day: day of measurement (0 or 7)

		\df_tidy_fluorimetry.csv
			mean_FvFm: mean Fv/Fm
			se_Fv/Fm: se Fv/Fm
			strain: strain identifier
			lifestyle: free-living (FL) vs. symbiotic (S)
			tech_rep: technical replicate

		\phylotypes.csv
			strain: strain identifier
			phylotype: clade based on phylogeny
			strain_long: assembly identifier from SSU phylogeny

		\brms_model_86.rds -- brms model for sugar release for all strains in SSU tree

		\brms_model_92.rds -- brms model for sugar release for all strains in plastid tree
		\brms_model.rds -- brms model for sugar release for all strains

		\mult_df_combo.csv -- data frame for multivariate analysis, means by strain derived from raw data frames for each phenotype


\phylogenetics -- Everything related to phylogenetic analysis of genomic data
(Note that intermediate files as well as raw genomic data have been omitted)

field_phylo_publish.Rproj -- RProject file for data analysis
labbook.md -- Notes on code for plastid tree

	\Data
		\Input
			\marker_trees -- treefiles for each plastid marker
			\blast_marker_summary.tsv -- summary of blast marker metrics
			\master_assembly_qc_table.tsv -- summary of plastid assembly metrics
			\phylogroups.csv -- phylogroups as judged from the final tree
			\plastid_tree4.treefile -- final treefile for tree included in publication
		\Output
			\best_samples.txt -- selected samples
			\core_markers.txt -- selected markers
			\final_sample_list.csv -- final samples list
		\Plots
			\tree_rotted_ponds_phylogroups.png
		\Scripts -- Scripts in order for completing phylogenetic analysis. See lab book.md for notes on what they were used for.
			\01-fastqc.slurm
			\02-multiqc.slurm
			\03-make_manifest.py
			\04-all_getorganelle.slurm
			\05-get_assembly_stats.py
			\06-add_metadata_to_assembly_stats.py
			\07-make_final_qc_assembly_stats_table.py
			\08-QCInspection.R
			\09-blast_many_markers.slurm
			\10-MarkerRecovery.R
			\11-blast_many_markers_references.slurm
			\12-extract_marker_fastas.slurm
			\13-add_more_extra_marker_fastas.slurm
			\14-align_markers.slurm
			\15-trim_markers.slurm
			\16-iqtree.slurm
			\17-marker_trees.slurm
			\18-Trees.R		
			
		\SSU_tree
			\final_alignments.mase -- alignments used to generate tree with SeaView for SSU tree by Tom Richards based on assemblies by Guy Leonard. The resulting tree was exported as pdf and manually edited in Illustrator.