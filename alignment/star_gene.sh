#!/bin/sh
#SBATCH -J star_align
#SBATCH --time=2:30:00
#SBATCH -p batch
#SBATCH -N 1 #nodes
#SBATCH -n 32 #tasks (can stand in for cores)
#SBATCH --mem=96gb #Memory requested
#SBATCH --output=sh_star_align.%j.%N.out
#SBATCH --error=sh_star_align.%j.%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=randall.ellis@tufts.edu

module load STAR

# Set variable for hg38 genomic index and output directory

INDEX=../human/star_index
OUTPUT_DIR=../output
WORKING_DIR=../fastq   # contains all of your fastqs

# Command line for STAR

cd ${WORKING_DIR}

#STAR assumes single-end if only one file used
for subject in *R1_001.fastq.gz
do
  SAMPLE=$(echo ${subject} | sed "s/_R1_001.fastq.gz//")
  echo ${SAMPLE}_R1_001.fastq.gz 
  STAR --genomeDir ${INDEX} \
  --runThreadN 32 \
  --readFilesIn ${SAMPLE}_R1_001.fastq.gz \
  --readFilesCommand zcat \
  --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE} \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts 
done
 
