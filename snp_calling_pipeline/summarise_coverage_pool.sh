#!/bin/bash
#SBATCH -J sumarise_poolseq_datasets
#SBATCH -o %x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=0-168:00:00
#SBATCH --account=def-yeaman
#SBATCH --mail-type=ALL
#SBATCH --mail-user=james.whiting@ucalgary.ca


# Loop over any...
datasets=(rellstab_Aalp rellstab_Ahal rellstab_Cres)
for dataset in "${datasets[@]}"
do

# Set dir
POOLSEQ_DIR=~/pool_snpcalling/$dataset/parentdir/

# Run the summarise coverage script
source $POOLSEQ_DIR/bash_variables
python ~/pipeline/98_get_read_stats.py $POOLSEQ_DIR 16

done
