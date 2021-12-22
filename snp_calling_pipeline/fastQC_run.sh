#!/bin/bash
#SBATCH -J fastQC_run
#SBATCH -o %J.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=0-12:00:00
#SBATCH --account=def-yeaman
#SBATCH --mail-type=ALL
#SBATCH --mail-user=james.whiting@ucalgary.ca

module load fastqc

MAIN=/home/jimw91/scratch/pool_snpcalling
DATASET=rellstab_Aalp

cd $MAIN/$DATASET/04_raw_data

parallel "fastqc {}" ::: $(ls *_1.fastq.gz)
