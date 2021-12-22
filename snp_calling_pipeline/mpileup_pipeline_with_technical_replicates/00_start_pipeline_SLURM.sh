#!/bin/bash
# Submit scripts from the $SPECIES_DIR directory

# Variables
MAIN=/home/jimw91/wgs_snpcalling/
DATASET=weigel_Athaliana_IBE
SPECIES_DIR=$MAIN/$DATASET
cd $SPECIES_DIR

# Point to scripts
PIPE_DIR=~/RepAdapt/snp_calling/general_scripts/snpcalling_pipeline_jw_multisamples/

# Set metadata
DATATABLE=$SPECIES_DIR/02_info_files/datatable.txt

# Set Email for slurm reports
EMAIL=james.whiting@ucalgary.ca

# Set ComputeCanada account
CC_ACCOUNT=def-yeaman

# How many samples are there?
FASTQ_N=$( ls $SPECIES_DIR/04_raw_data/*fastq.gz | wc -l )
FILE_ARRAY=$(( $FASTQ_N / 2 ))

'''
##########################
# Part 1 of the pipeline #
##########################
'''

# Trim
job01=$(sbatch --account=$CC_ACCOUNT \
    --array=1-${FILE_ARRAY} \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --parsable \
    $PIPE_DIR/01_fastp.sh)

# Index reference & Align reads to reference
# Note - If fastq files include .1 and .2 suffixes, bwa will fail. Lines in script 02 can be commented out to handle this
job02=$(sbatch --account=$CC_ACCOUNT  \
    --array=1-${FILE_ARRAY} \
    --dependency=afterok:$job01 \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --parsable \
    $PIPE_DIR/02_bwa_alignments.sh)

# Collect sample data metrics
job03=$(sbatch --account=$CC_ACCOUNT  \
    --array=1-${FILE_ARRAY} \
    --dependency=afterok:$job02 \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --parsable \
    $PIPE_DIR/03_collect_metrics.sh)

'''
##########################
# Part 2 of the pipeline #
##########################
'''

# Remove duplicates
job04=$(sbatch --account=$CC_ACCOUNT  \
    --array=1-${FILE_ARRAY} \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --parsable \
    $PIPE_DIR/04_remove_duplicates.sh)


# Change bam files RG
job05=$(sbatch --account=$CC_ACCOUNT  \
    --array=1-${FILE_ARRAY} \
    --dependency=afterok:$job04 \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --parsable \
    $PIPE_DIR/05_change_RG.sh)

'''
##########################
# Part 3 of the pipeline #
##########################
'''

# New step here now to merge BAM files over samples...
# Count the number of unique samples in $DATATABLE
SAMPLE_ARRAY=$(cut -f14 $DATATABLE | sort | uniq | wc -l)
job05b=$(sbatch --account=$CC_ACCOUNT  \
    --array=1-${SAMPLE_ARRAY} \
    --dependency=afterok:$job05 \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --parsable \
    $PIPE_DIR/05b_merge_bams.sh)

# Realign around indels...
job06=$(sbatch --account=$CC_ACCOUNT  \
    --array=1-${SAMPLE_ARRAY} \
    --dependency=afterok:$job05b \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --parsable \
    $PIPE_DIR/06_gatk_realignments.sh)

# Stats of final bam files...
job06b=$(sbatch --account=$CC_ACCOUNT  \
    --array=1-${SAMPLE_ARRAY} \
    --dependency=afterok:$job06 \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --parsable \
    $PIPE_DIR/06b_collect_final_metrics.sh)

'''
##########################
# Part 4 of the pipeline #
##########################
'''

##### Set up scaffold input files to fit with Compute Canada max jobs
# How many scaffolds are in the genome...
SCAFF_N=$(cat $SPECIES_DIR/03_genome/*fai | wc -l)
SPLIT_N=200

# Split these over $SPLIT_N jobs...
cut -f1 $SPECIES_DIR/03_genome/*fai | shuf > 02_info_files/all_scafs.txt
if [[ $SCAFF_N -gt $SPLIT_N ]]
then
  split -l$((`wc -l < 02_info_files/all_scafs.txt`/${SPLIT_N})) 02_info_files/all_scafs.txt 02_info_files/all_scafs.split. -da 4 --additional-suffix=".pos"
else
  split -l$((`wc -l < 02_info_files/all_scafs.txt`/${SCAFF_N})) 02_info_files/all_scafs.txt 02_info_files/all_scafs.split. -da 4 --additional-suffix=".pos"
fi

# Set SNP-calling array over these scaffold clusters...
SCAFF_ARRAY=$(ls 02_info_files/all_scafs*pos | wc -l)

##########################
##### Make some metadata
#### Bam List...
ls 06_bam_files/*realigned.bam > 02_info_files/bammap.txt

#### Ploidy file...
# Ploidy information is pulled from the second column of datatable.txt and should be available with SRA/ENA metadata.
#awk 'BEGIN {OFS="\t"} {print $7,$2}' 02_info_files/datatable.txt | sort -nk1 | uniq > 02_info_files/ploidymap.txt
# Ploidy information is built from bam list
rm -f 02_info_files/ploidymap.txt
for BAM in $(cat 02_info_files/bammap.txt)
do
  BAM_CLEAN=$(basename $BAM | sed 's/.realigned.bam//g')
  BAM_PLOIDY=$(grep -w $BAM_CLEAN 02_info_files/datatable.txt | cut -f4 | head -n1)
  echo -e "$BAM_CLEAN\t$BAM_PLOIDY" >> 02_info_files/ploidymap.txt
done

##########################
# Call SNPs - Mpileup runs first and filtering starts based on the dependency
export DATASET=$DATASET
job07=$(sbatch --account=$CC_ACCOUNT \
    --array=1-${SCAFF_ARRAY} \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --export DATASET \
    --parsable \
    $PIPE_DIR/07_mpileup.sh)

# # Filter the SNPs
# export DATASET=$DATASET
# job08=$(sbatch --account=$CC_ACCOUNT \
#     --dependency=afterok:$job07 \
#     --array=1-${SCAFF_ARRAY} \
#     --mail-type=ALL \
#     --mail-user=$EMAIL \
#     --export DATASET \
#     --parsable \
#     $PIPE_DIR/08_scaffoldVCF_filtering.sh)

# Filter the SNPs
export DATASET=$DATASET
job08=$(sbatch --account=$CC_ACCOUNT \
    --array=1-${SCAFF_ARRAY} \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --export DATASET \
    --parsable \
    $PIPE_DIR/08b_relaxed_scaffoldVCF_filtering.sh)

# Concatenate the per-scaffold VCFs to a single VCF
export DATASET=$DATASET
sbatch --account=$CC_ACCOUNT \
    --mail-type=ALL \
    --dependency=afterok:$job08 \
    --mail-user=$EMAIL \
    --export DATASET \
    $PIPE_DIR/09_concat_VCFs.sh

##########################
