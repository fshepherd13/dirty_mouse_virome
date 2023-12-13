#!/bin/bash
#SBATCH -J "Create snakemake 7.20.0 env"
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=1g
#SBATCH --time=2:00:00
#SBATCH -M agate
#SBATCH -p langlois-node1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sheph085@umn.edu

module load mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake=7.26.0
