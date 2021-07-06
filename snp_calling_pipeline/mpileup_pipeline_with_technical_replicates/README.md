# RepAdapt mpileup SNP Calling Pipeline for WGS/Capture
*This version of the pipeline allows for technical replicates to be handled on the basis of a metadata file, with sample IDs provided and SRA/ENA codes pointing to the same sample for merging at step 5b.*

## Setup
### Directory Structure
The pipeline assumes that all individuals are from the same species and assumes the following directory and subdirectory structure that can be set up with:
```
# Make a parent directory for species
SPECIES_DIR=name_of_your_species
mkdir $SPECIES_DIR

# Make subdirectories
cd $SPECIES_DIR
mkdir 00_archive  02_info_files  04_raw_data	06_bam_files  08_filtered_VCFs  99_metrics	01_scripts  03_genome	 05_trimmed_data  07_raw_VCFs 98_log_files
```

### Raw Data
Throughout here I've assumed starting from raw SRA/ENA read files, for e.g. `SRR10335077_1.fastq.gz`. Will rename individuals to meaningful sample names in the final VCF, so keep a metadata linking accessions to samples. To identify readpairs then, usually these are listed as `ls *1.fastq.gz`. If for whatever reason this doesn't match the format of your raw read data, then those parts of the scripts that use `ls *1.fastq.gz` will need editing.

All raw data should be placed in `04_raw_data/`

If starting from BAM files, these can be stored in `06_bam_files/`

### Genome
Also need to prepare the genome files and move into `03_genome`
```
# Make sure genome is suffixed with one of "".fasta .fa .fasta.gz or .fa.gz" as this is used later to identify the genome file within the subdir
GENOME=your_genome_file.fasta

# bwa index
bwa index $GENOME

# samtools index
samtools faidx $GENOME

# gatk dictionary with picard
java -jar $EBROOTPICARD/picard.jar CreateSequenceDictionary R=$GENOME O=${GENOME%.*}.dict
```
All genome files can either be moved into the `03_genome` subdirectory, or if using the same genome for several species, symlinking them may be simpler, ie. `ln -S /path/to/genome* $SPECIES_DIR/03_genome/`

### DataTable
This is an important metadata file that should include the following columns (even if empty) and must be formatted as such:

|sample_name|library_name|pool_name|ploidy|file_name_r1          |file_name_r2          |adaptor_1|adaptor_2|ref|rgid                                      |rglb|rgpl    |rgpu|rgsm                          |
|-----------|------------|---------|------|----------------------|----------------------|---------|---------|---|------------------------------------------|----|--------|----|------------------------------|
|SRR10340174|            |         |2     |SRR10340174_1.fastq.gz|SRR10340174_2.fastq.gz|         |         |   |Eucalyptus_albens_Nangar_J361b_SRR10340174|    |ILLUMINA|    |Eucalyptus_albens_Nangar_J361b|
|SRR10340175|            |         |2     |SRR10340175_1.fastq.gz|SRR10340175_2.fastq.gz|         |         |   |Eucalyptus_albens_Nangar_J361a_SRR10340175|    |ILLUMINA|    |Eucalyptus_albens_Nangar_J361a|
|SRR10340176|            |         |2     |SRR10340176_1.fastq.gz|SRR10340176_2.fastq.gz|         |         |   |Eucalyptus_albens_Nangar_J361_SRR10340176 |    |ILLUMINA|    |Eucalyptus_albens_Nangar_J361 |
|SRR10340177|            |         |2     |SRR10340177_1.fastq.gz|SRR10340177_2.fastq.gz|         |         |   |Eucalyptus_albens_Nangar_J361_SRR10340177 |    |ILLUMINA|    |Eucalyptus_albens_Nangar_J361 |

### Ploidy
Ploidy information should be added to the metadata file, from which the ploidymap is derived automatically.

---

## Running
The script called `00_start_pipeline_SLURM.sh` includes `sbatch` commands for running each of the steps of the pipeline. This script should be used interactively, so to run each step, first set the general variables, for example:
```
# Variables
MAIN=/home/jimw91/RepAdapt/snp_calling
DATASET=murray_Emol
SPECIES_DIR=$MAIN/$DATASET
cd $SPECIES_DIR

# Point to scripts
PIPE_DIR=$MAIN/general_scripts/snpcalling_pipeline_jw

# Set Email for slurm reports
EMAIL=james.whiting@ucalgary.ca

# Set ComputeCanada account
CC_ACCOUNT=def-yeaman

# How many samples are there?
FASTQ_N=$( ls $SPECIES_DIR/04_raw_data | wc -l )
FILE_ARRAY=$(( $FASTQ_N / 2 ))
```

Each step of the pipeline can then just be run as batch array jobs using the above variables to receive status emails and ensure that each step of the pipeline is starting from the correct working directory.

### Dependencies
I have also added job dependencies so that scripts are run in groups. These are currently grouped as:
Submit scripts 1 + 2 + 3
Submit scripts 4 + 5 + 6
Submit scripts 7 + 8 + 9

These groupings just reflect points where it's useful to check outputs, and also if running lots of individuals/lots of scaffolds compute canada limits to 1000 jobs which prevents running them all as a single group. There's also some metadata to make in between scripts 6 and 7, so splitting into scripts 1-6 and then 7-9 may also be sensible.

To remove dependencies, just remove lines that include `--dependency=`.

---

## Debugging
Log files are written to `98_log_files/`, and are named according to the format of ${JOB_NAME}_${ARRAY_NUMBER}_array${SLURM_ARRAY_TASK_ID}. So to trace back problems, identify the problematic file in the array order (will match `ls` order), and check relevant log file.

---

## Notes
For script 2 (bwa), depending on how fastq files are downloaded from SRA the fastq files may have added a .1 or .2 suffix to the read pair names. This causes bwa to fail as it can't match read pairs from R1 and R2 fastqs. The following lines in script 2 can be uncommented and used to fix this by just removing these suffixes. By default these are commented out:
```
#######################################
# Uncomment these lines in rare cases where SRA download has suffixed R1 and R2 fastq headers with .1 and .2, which errors out bwa. This removes them
###
mv $RAWDATAFOLDER/$file1 $RAWDATAFOLDER/${name}.R1.trimmed_dirty.fastq.gz
mv $RAWDATAFOLDER/$file2 $RAWDATAFOLDER/${name}.R2.trimmed_dirty.fastq.gz
zcat $RAWDATAFOLDER/${name}.R1.trimmed_dirty.fastq.gz | sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" | gzip > $RAWDATAFOLDER/$file1
zcat $RAWDATAFOLDER/${name}.R2.trimmed_dirty.fastq.gz | sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" | gzip > $RAWDATAFOLDER/$file2
#######################################
```
