# RepAdapt WGS/Capture SNP-Calling pipeline
This directory includes the scripts required to call SNPs from raw fastq files.

The scripts are numbered 01-09, with a wrapper for running jobs on a SLURM-based HPC.

Each script is described more below. 

Steps 1-6 are run over individual files as batch arrays.
Steps 7-8 are run over groups of scaffolds of size N as batch arrays.
Step 9 is run as a single job.

## Metadata and Genome preparation
Each dataset requires a metadata file that must be in the following format. It must be stored as `02_info_files/datatable.txt`

| #SRA       |   |   | ploidy | R1_file                | R2_file                |   |   |   | RG                           |   | instrument |   | individual        | r1_ftp                                                                   | r2_ftp                                                                   | r1_md5                           | r2_md5                           |
|------------|---|---|--------|------------------------|------------------------|---|---|---|------------------------------|---|------------|---|-------------------|--------------------------------------------------------------------------|--------------------------------------------------------------------------|----------------------------------|----------------------------------|
| SRR2744682 |   |   | 2      | SRR2744682_R1.fastq.gz | SRR2744682_R2.fastq.gz |   |   |   | Ptremula_SwAsp001_SRR2744682 |   | ILLUMINA   |   | Ptremula_SwAsp001 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR274/002/SRR2744682/SRR2744682_1.fastq.gz | ftp.sra.ebi.ac.uk/vol1/fastq/SRR274/002/SRR2744682/SRR2744682_2.fastq.gz | a329816ed2d49ca6d05bf3efc5801024 | f074544db22d60f0ca8ea573f6ac7bdf |
| SRR4301375 |   |   | 2      | SRR4301375_R1.fastq.gz | SRR4301375_R2.fastq.gz |   |   |   | Ptremula_SwAsp003_SRR4301375 |   | ILLUMINA   |   | Ptremula_SwAsp003 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR430/005/SRR4301375/SRR4301375_1.fastq.gz | ftp.sra.ebi.ac.uk/vol1/fastq/SRR430/005/SRR4301375/SRR4301375_2.fastq.gz | eb2a9b1b1e2a88542342dd3b6f3445a3 | 2d18328dea665ab74d2a77e77799372d |
| SRR4301376 |   |   | 2      | SRR4301376_R1.fastq.gz | SRR4301376_R2.fastq.gz |   |   |   | Ptremula_SwAsp004_SRR4301376 |   | ILLUMINA   |   | Ptremula_SwAsp004 | ftp.sra.ebi.ac.uk/vol1/fastq/SRR430/006/SRR4301376/SRR4301376_1.fastq.gz | ftp.sra.ebi.ac.uk/vol1/fastq/SRR430/006/SRR4301376/SRR4301376_2.fastq.gz | 567a4448a76217102fdecc1b78020aec | 20f3166ac143267c9ea987f71b6c811c |

The reference genome must also be pre-processed by indexing with bwa and samtools and preparing the dictionary files for gatk.
This can be done using the command: `sh prepare_genome.sh $ref_genome`

The pipeline also expects a directory structure that can be made using the following command
```
mkdir 00_archive  01_scripts  02_info_files  03_genome  04_raw_data  05_trimmed_data  06_bam_files  07_raw_VCFs  08_filtered_VCFs  09_final_vcf  98_log_files  99_metrics  99_metrics_merged
```

## General Scripts
A few general scripts are included which perform the following tasks:
* `ena_wget_downloader.sh` - Uses the metadata file to download raw data from ena and produce and check md5
* `fastq_dump_parallel.sh` - An alternative approach to downloading raw fastq data using the sra-toolkit (tended to be a bit buggy)
* `concat_bam_stats.sh` - Point at a directory which includes per-bam summary stats to concatenate to a single file
* `prepare_genome.sh` - Point at a reference genome to perform necessary indexing and additional genome files

## Step 0 - Wrapper for running Steps 1-9 on an HPC
This script allows the user to set metadata regarding where pipeline scripts and the working directory is.
Also set email for notifications and general HPC info
Each script is then submitted as an individual job with dependencies on the prior job completing successfully before running.

## Step 1 - Read Trimming
Read trimming is performed using fastp.

* Script: `01_fastp.sh`
* Inputs: Raw fastq files
* Outputs: Trimmed fastq files

## Step 2 - Alignment to Genome
Alignment is performed using bwa-mem.
Also includes a QC on trimmed fastq to make sure read headers are not malformed, and attempts to resolve if they are.

* Script: `02_bwa_alignments.sh`
* Inputs: Trimmed fastq files
* Outputs: Sorted and indexed BAM files.

## Step 3 - First run of BAM quality metrics
BAM quality metrics calculcated using Picard:
(`CollectAlignmentSummaryMetrics`, `CollectInsertSizeMetrics`,`CollectWgsMetricsWithNonZeroCoverage`)

* Script: `03_collect_metrics.sh`
* Inputs: Sorted/Indexed BAMs
* Outputs: Alignment, insert size and WGS metrics, outputted to `99_metrics/`

## Step 4 - Remove duplicates
Duplicates are marked and removed using Picard's `MarkDuplicates`

* Script: `04_remove_duplicates.sh`
* Inputs: Sorted/Indexed BAMs
* Outputs: Deduplicated BAMs (`*dedup.bam`)

## Step 5 - Change Read Groups and Merge Replicates
Modifies read groups and then merges over technical replicates of the same individual
Read groups are added using Picard's `AddOrReplaceReadGroups` and merging is done with `MergeSamFiles`

Note: Even if there aren't any technical replicates, the merging step should still be run to maintain filenames.

* Scripts: `05_change_RG.sh` and `05b_merge_bams.sh`
* Inputs: Deduplicated BAMs (`*.dedup.bam`)
* Outputs: Merged BAMs (`*.merged.bam`)

## Step 6 - Realignment around indels and final BAM quality checks
Performs indel realignment by first generating interval files using GATK's `RealignerTargetCreator`.
Indel realignment is then performed with `IndelRealigner`.
Final quality metrics are then calculated over the final BAM files as in step 3.

Note: These software are only available in versions of GATK prior to v4.

* Scripts: `06_gatk_realignments.sh` and `06b_collect_final_metrics.sh`
* Inputs: Merged BAMs (`*.merged.bam`)
* Outputs: Realigned BAMs (`*.realigned.bam`)

Note: These are deposited in `99_metrics_merged/`. 
Can then run `concat_bam_stats.sh 99_metrics_merged/ DATASET_NAME` to summarise to a single text file.

## Step 7 - SNP calling through mpileup
Performs SNP-calling using genotype-likelihoods.
The genome is split into groups of scaffolds, and these groups are run in parallel through BCFtools' `mpileup` and `call`.

* Script: `07_mpileup.sh`
* Inputs: All realigned BAMs
* Outputs: Per-scaffold unfiltered VCFs

## Step 8 - Scaffold VCF filtering
Each scaffold is filtered according to a standardised set of quality filters
Filtering is performed using VCFtools, specifically the following:
```
--minQ 30 \
--minGQ 20 \
--minDP 5 \
--max-alleles 2 \
--max-missing 0.7 \
```

* Script: `08_scaffoldVCF_filtering.sh`
* Inputs: Per-scaffold unfiltered VCFs
* Outputs: Per-scaffold quality-filtered VCFs

## Step 9 - Final VCF filtering and concatenation
Script concatenates all per-scaffold filtered VCFs and performs some final maf-filtering using BCFtools

* Script: `09_concat_VCFs.sh`
* Inputs: Per-scaffold quality-filtered VCFs
* Outputs: Three final VCFs with varying maf-filtering of: None, 1% and 5%:
```
*_full_concatened.vcf
*_full_concatened_maf01.vcf.gz
*_full_concatened_maf05.vcf.gz
```