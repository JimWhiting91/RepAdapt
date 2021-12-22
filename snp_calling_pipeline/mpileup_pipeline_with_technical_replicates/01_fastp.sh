#!/bin/bash
# 6 CPUs
# 10 Go

#SBATCH -J 01.fastp
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
##SBATCH --partition=cpu2019,cpu2013,lattice,apophis-bf
#SBATCH --mem=10G
#SBATCH --time=0-06:00:00

# Load up fastp
module load StdEnv/2020 fastp/0.20.1

#cd $SLURM_SUBMIT_DIR

##Keep some info. about the run/script
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)

# Variables
INDIR="04_raw_data"
OUTDIR="05_trimmed_data"
LOG="98_log_files"
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
mkdir $OUTDIR/01_reports

# Make a log file to the species log directory
cp $SCRIPT $LOG/"$TIMESTAMP"_"$NAME"_%J

# Pull file from the FASTP_ARRAY
input_file=$(cut -f1 02_info_files/datatable.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")

# Run over file
    #input_file=$(echo "$file" | perl -pe 's/_R1.*\.fastq.gz//')
    output_file=$(basename "$input_file")
    echo "Still working for you... Cleaning: $input_file"

    fastp -w ${SLURM_CPUS_PER_TASK} \
        -i $INDIR/${input_file}_R1.fastq.gz \
        -I $INDIR/${input_file}_R2.fastq.gz \
        -o $OUTDIR/"$output_file".R1.trimmed.fastq.gz \
        -O $OUTDIR/"$output_file".R2.trimmed.fastq.gz \
        -j $OUTDIR/01_reports/"$output_file".json \
        -h $OUTDIR/01_reports/"$output_file".html \
        &> "$LOG"/01_fastp_"$output_file"_"$TIMESTAMP".out
