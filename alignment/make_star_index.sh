#!/bin/sh
#SBATCH -J make_star_index
#SBATCH --time=04:00:00
#SBATCH -p batch
#SBATCH -N 1 #nodes
#SBATCH -n 16 #tasks (can stand in for cores)
#SBATCH --mem=96gb #Memory requested
#SBATCH --output=sh_make_star_index.%j.%N.out
#SBATCH --error=sh_make_star_index.%j.%N.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=randall.ellis@tufts.edu

module load STAR

STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir /cluster/tufts/levinlab/rellis01/ngs/genomes/human/star_index \
     --genomeFastaFiles /cluster/tufts/levinlab/rellis01/ngs/genomes/human/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
     --sjdbGTFfile /cluster/tufts/levinlab/rellis01/ngs/genomes/human/Homo_sapiens.GRCh38.110.gtf 
