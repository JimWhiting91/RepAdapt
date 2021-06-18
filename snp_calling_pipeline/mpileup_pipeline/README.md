# RepAdapt mpileup SNP Calling Pipeline for WGS/Capture

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

### Ploidy
If ploidy is known, then a ploidy file can be made and placed in `02_info_files/ploidymap.txt`. By default (because I don't have any ploidy info on the eucalypt samples), a basic diploid map is made assuming all samples are diploid:
```
ls 04_raw_data/*1.fastq.gz | xargs -n 1 basename | sed 's/_1.fastq.gz//g' > 02_info_files/ploidymap.txt
for ind in $(cut -f1 02_info_files/ploidymap.txt)
do
  sed -i "s/$ind/${ind}\t2/g" 02_info_files/ploidymap.txt
done
```

Here, column 1 is sample ID (e.g. SRR10335077) and column 2 is ploidy count.

If you're using known ploidy information, comment these out or just don't run them.

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
