#!/bin/bash
#SBATCH -J prepare_genome
#SBATCH -o %x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=0-12:00:00
#SBATCH --account=def-yeaman
#SBATCH -D .
#SBATCH --mail-type=ALL
#SBATCH --mail-user=james.whiting@ucalgary.ca

module load bwa samtools picard

GENOME=$1

cd $(dirname $GENOME)

# bwa index
bwa index $GENOME

# samtools index
samtools faidx $GENOME

# gatk dictionary with picard
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=$GENOME O=${GENOME%.*}.dict

