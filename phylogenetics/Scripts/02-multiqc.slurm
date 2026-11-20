#!/bin/bash --login

#SBATCH -p serial
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH --job-name=multiqc
#SBATCH -o multiqc_%j.out

# A script to summarise the fastqc results

module purge
module load apps/python/multiqc/1.12

PROJECT=PROJECT_FOLDER

mkdir -p $PROJECT/qc/multiqc

multiqc \
	$PROJECT/qc/fastqc_results \
	-o $PROJECT/qc/multiqc 
