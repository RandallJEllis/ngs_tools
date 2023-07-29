#BSUB -J STAR_align_paired_rsem  ## name whatever you want job to show up as
#BSUB -P acc_DADisorders
#BSUB -q premium
#BSUB -n 20
#BSUB -R span[hosts=1]
#BSUB -R rusage[mem=16000]
#BSUB -W 5:00     ## For most RNA-seq, probably max 1 hour needed per subject, sometimes asking for less $
#BSUB -o STAR_align_paired_rsem.stdout  ## output files -- name them accordingly
#BSUB -eo STAR_align_paired_rsem.stderr  ## output files --name accordingly
#BSUB -L /bin/bash

ml star

# Set variable for hg38 genomic index and output directory

INDEX=../genomes/rat/starIndex/
OUTPUT_DIR=../output/transcriptome/TRANSCRIPTOME_FASTQ_2023
WORKING_DIR=../FASTQ_GENEWIZ_2023   # contains all of your fastqs

# Command line for STAR

cd ${WORKING_DIR}

#STAR assumes single-end if only one file used
for subject in *R1_001.fastq
do
  SAMPLE=$(echo ${subject} | sed "s/_R1_001.fastq//")
  echo ${SAMPLE}_R1_001.fastq ${SAMPLE}_R2_001.fastq
  STAR --genomeDir ${INDEX} \
  --runThreadN 20 \
  --readFilesIn ${SAMPLE}_R1_001.fastq ${SAMPLE}_R2_001.fastq \
  --outFileNamePrefix ${OUTPUT_DIR}/${SAMPLE} \
  --outSAMtype BAM SortedByCoordinate \
  --quantMode TranscriptomeSAM
done
 
