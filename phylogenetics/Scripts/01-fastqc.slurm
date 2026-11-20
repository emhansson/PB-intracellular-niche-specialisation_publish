#!/bin/bash --login

#SBATCH -p multicore 
#SBATCH -n 4
#SBATCH -t 8:00:00
#SBATCH --job-name=fastqc
#SBATCH -o fastqc_%j.out

# A script to run fastqc
# Note that it also creates the fastqc_temp directory with symlinks with unique names
# for the files. This is also used for the downstream organelle assembly because the directory
# structure/naming for the raw files is ridiculous.

module purge
module load apps/binapps/fastqc/0.12.1

PROJECT=PROJECT_FOLDER

mkdir -p $PROJECT/qc/fastqc_results
mkdir -p $PROJECT/qc/fastqc_temp

find $PROJECT/ch_like \
\( -name "*_R1_001.fastq.gz" -o -name "*_R2_001.fastq.gz" \) \
| while read file; do

	sampledir=$(basename "$(dirname "$file")")
	filename=$(basename "$file")

	unique="${sampledir}_${filename}"

	ln -s "$file" "$PROJECT/qc/fastqc_temp/$unique"

done

find $PROJECT/qc/fastqc_temp \
	-name "*.fastq.gz" \
	| xargs -n 1 -P 4 fastqc -t 1 -o $PROJECT/qc/fastqc_results
