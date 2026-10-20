#!/bin/bash

#SBATCH --job-name=vcf_stats

#SBATCH --partition=general

#SBATCH --nodes=1

#SBATCH --cpus-per-task=1

#SBATCH --mem=5G

#SBATCH --time=1:00:00

#SBATCH --output=./logs/%x.%j.out

#SBATCH --error=./logs/%x.%j.err

#SBATCH --mail-type=ALL

#SBATCH --mail-user=wreganlo@uvm.edu

cd ~/projects/eco_genomics_2026/population_genetics/myscripts
module load gcc
module load vcftools



vcftools --gzvcf Lottia_gigantea_biallelic.vcf.gz --missing-indv --out missing_data

vcftools --gzvcf Lottia_gigantea_biallelic.vcf.gz --depth --out depth_data
