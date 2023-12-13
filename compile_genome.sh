#!/bin/bash
#SBATCH -J "Generate reference genome"
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --mem=50g
#SBATCH --time=3:00:00
#SBATCH -M mesabi
#SBATCH -p amdsmall
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sheph085@umn.edu

module load star/2.7.1a
#Define directory where you want your reference genome to go
DIR="/home/langlois/shared/ref_genomes/star_indices/mus_musculus"

mkdir -p $DIR
cd $DIR

# Define locations online (Ensembl or other) for the reference genome fasta file and GTF annotation file
GENOME="https://ftp.ensembl.org/pub/release-109/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"
GTF="https://ftp.ensembl.org/pub/release-109/gtf/mus_musculus/Mus_musculus.GRCm39.109.gtf.gz"

#Download with wget and decompress
wget -O - $GENOME | gunzip -c > genome.fa
wget -O - $GTF | gunzip -c > annotation.gtf

#Create file that logs which assembly and annotation was downloaded for future reference
echo "Downloaded $GENOME as the fasta file and $GTF as the annotation file." > genome_info.txt

#run STAR genomeGenerate
STAR --runThreadN 16 \
--runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles genome.fa \
--sjdbGTFfile annotation.gtf \
--sjdbOverhang 149