#!/bin/bash
# 1 CPU
# 30 Go

#SBATCH -J 04.Duplicates
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem 200G
#SBATCH --time=00-01:00:00

# Load required modules
module load picard java

# Global variables
PICARD=$EBROOTPICARD/picard.jar
MARKDUPS="MarkDuplicates"
ALIGNEDFOLDER="06_bam_files"
METRICSFOLDER="99_metrics"

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"

export JAVA_TOOL_OPTIONS="-Xms2g -Xmx50g "
export _JAVA_OPTIONS="-Xms2g -Xmx50g "

# Fetch filename from the array
sample_name=$(cut -f1 02_info_files/datatable.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")
file=${sample_name}.sorted.bam

    echo "DEduplicatING sample $file"

    java -jar $PICARD $MARKDUPS \
        INPUT=$ALIGNEDFOLDER/$file \
        OUTPUT=$ALIGNEDFOLDER/${sample_name}.dedup.bam \
        METRICS_FILE=$METRICSFOLDER/${sample_name}_DUP_metrics.txt \
        VALIDATION_STRINGENCY=SILENT \
        REMOVE_DUPLICATES=true

echo " >>> Cleaning a bit...
"
echo "DONE! Go check your files."
