#!/bin/bash
# Submit scripts from the $SPECIES_DIR directory

# Variables
MAIN=/home/jimw91/RepAdapt/snp_calling
DATASET=murray_Emol
SPECIES_DIR=$MAIN/$DATASET
cd $SPECIES_DIR

# Point to scripts
PIPE_DIR=$MAIN/general_scripts/snpcalling_pipeline_jw

# Set Email for slurm reports
EMAIL=james.whiting@ucalgary.ca

# How many samples are there?
FASTQ_N=$( ls $SPECIES_DIR/04_raw_data | wc -l )
FILE_ARRAY=$(( $FASTQ_N / 2 ))

'''
##########################
# Part 1 of the pipeline #
##########################
'''

# Trim
sbatch --account=def-yeaman \
    --array=1-${FILE_ARRAY} \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    $PIPE_DIR/01_fastp.sh

# Index reference & Align reads to reference
# Note - If fastq files include .1 and .2 suffixes, bwa will fail. Lines in script 02 can be commented out to handle this
sbatch --account=def-yeaman  \
    --array=1-${FILE_ARRAY} \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    $PIPE_DIR/02_bwa_alignments.sh

'''
##########################
# Part 2 of the pipeline #
##########################
'''

# Collect sample data metrics
sbatch --account=def-yeaman  \
    --array=1-${FILE_ARRAY} \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    $PIPE_DIR/03_collect_metrics.sh

# Remove duplicates
sbatch --account=def-yeaman  \
    --array=1-${FILE_ARRAY} \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    $PIPE_DIR/04_remove_duplicates.sh

'''
##########################
# Part 3 of the pipeline #
##########################
'''

# Change bam files RG
sbatch --account=def-yeaman  \
    --array=1-${FILE_ARRAY} \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    $PIPE_DIR/05_change_RG.sh

# Realign around indels...
sbatch --account=def-yeaman  \
    --array=1-${FILE_ARRAY} \
    -D $SPECIES_DIR \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    $PIPE_DIR/06_gatk_realignments.sh

'''
##########################
# Part 4 of the pipeline #
##########################
'''

##### Set up scaffold input files to fit with Compute Canada max jobs
# How many scaffolds are in the genome...
SCAFF_N=$(cat $SPECIES_DIR/03_genome/*fai | wc -l)
SPLIT_N=200

# Split these over 500 jobs...
cut -f1 $SPECIES_DIR/03_genome/*fai > 02_info_files/all_scafs.txt
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
# BY DEFAULT, PLOIDY HERE IS SET AS DIPLOID. IF THIS IS NOT THE CASE, EDIT IN SCRIPT.
# If ploidy information is available per individual, save a different version of 02_info_files/ploidymap.txt
ls 04_raw_data/*1.fastq.gz | xargs -n 1 basename | sed 's/_1.fastq.gz//g' > 02_info_files/ploidymap.txt
for ind in $(cut -f1 02_info_files/ploidymap.txt)
do
  sed -i "s/$ind/${ind}\t2/g" 02_info_files/ploidymap.txt
done

##########################
# Call SNPs - Mpileup runs first and filtering starts based on the dependency
export DATASET=$DATASET
job07=$(sbatch --account=def-yeaman \
    --array=1-${SCAFF_ARRAY} \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --export DATASET \
    --parsable \
    $PIPE_DIR/07_mpileup.sh)

# Filter the SNPs
export DATASET=$DATASET
sbatch --account=def-yeaman \
    --dependency=afterok:$job07 \
    --array=1-${SCAFF_ARRAY} \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --export DATASET \
    $PIPE_DIR/08_scaffoldVCF_filtering.sh

# Concatenate the per-scaffold VCFs to a single VCF
export DATASET=$DATASET
sbatch --account=def-yeaman \
    --mail-type=ALL \
    --mail-user=$EMAIL \
    --export DATASET \
    $PIPE_DIR/09_concat_VCFs.sh

##########################
