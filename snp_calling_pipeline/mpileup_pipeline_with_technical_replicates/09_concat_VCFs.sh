#!/bin/bash

#SBATCH --job-name="09.concatVCF"
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00-00:30:00

module load vcftools
module load StdEnv/2020 intel/2020.1.217 tabix/0.2.6

cd $SLURM_SUBMIT_DIR

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"
echo $SCRIPT
cp $SCRIPT ${LOG_FOLDER}/${TIMESTAMP}_${NAME}

# Variables
VCF="07_raw_VCFs"
FILTVCF="08_filtered_VCFs"

begin=`date +%s`

# Concatenate all the scaffold-VCF files into one global VCF file
vcf-concat $(ls -1 $FILTVCF/*_filtered.vcf.gz | perl -pe 's/\n/ /g') > ${FILTVCF}/${DATASET}_full_concatened.vcf && bgzip ${FILTVCF}/${DATASET}_full_concatened.vcf
# 2> "$LOG_FOLDER"/09_concatVCF_"$TIMESTAMP".log

echo "
DONE! Check you files"

end=`date +%s`
elapsed=`expr $end - $begin`
echo Time taken: $elapsed s
