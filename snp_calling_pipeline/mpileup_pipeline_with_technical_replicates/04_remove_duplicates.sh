#!/bin/bash
# 1 CPU
# 30 Go

#SBATCH -J 04.Duplicates
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 30G
#SBATCH --time=00-01:00:00

# Load required modules
module load picard java

# Global variables
PICARD=$EBROOTPICARD/picard.jar
MARKDUPS="MarkDuplicates"
ALIGNEDFOLDER="06_bam_files"
METRICSFOLDER="99_metrics"
# BAM_FILE="02_info_files/bam_split/bam0000"

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"

# Remove duplicates from bam alignments
# cat "$BAM_FILE" |
# while read file
# do

# Fetch filename from the array
file=$(ls $ALIGNEDFOLDER/*.sorted.bam | sed "${SLURM_ARRAY_TASK_ID}q;d" | xargs -n 1 basename)
sample_name=${file%.*.*}

    echo "DEduplicatING sample $file"

    java -jar $PICARD $MARKDUPS \
        INPUT=$ALIGNEDFOLDER/$file \
        OUTPUT=$ALIGNEDFOLDER/${sample_name}.dedup.bam \
        METRICS_FILE=$METRICSFOLDER/${sample_name}_DUP_metrics.txt \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true
# done 2> "$LOG_FOLDER"/04_duplicates_"$TIMESTAMP".log

echo " >>> Cleaning a bit...
"
#rm "$ALIGNEDFOLDER"/"$file"
echo "DONE! Go check your files."
