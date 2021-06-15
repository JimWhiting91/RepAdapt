#!/bin/bash
# 1 CPU
# 10 Go

#SBATCH -J 03.Metrics
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem 10G
#SBATCH --time=0-12:00:00

# Load modules
module load picard java

# Global variables
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/*{fasta,fa,fasta.gz,fa.gz} | xargs -n 1 basename)
ALIGNEDFOLDER="06_bam_files"
METRICSFOLDER="99_metrics"
PICARD=$EBROOTPICARD/picard.jar
ALIGN="CollectAlignmentSummaryMetrics"
INSERT="CollectInsertSizeMetrics"
COVERAGE="CollectWgsMetricsWithNonZeroCoverage"
# BAM_FILE="02_info_files/bam_split/bam0000"

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
SCRIPTNAME=$(basename $0)
LOG_FOLDER="98_log_files"
cp $SCRIPT $LOG_FOLDER/${TIMESTAMP}_${SCRIPTNAME}


# Load needed modules
# module load java
#module load samtools/1.8

# Run Picard Tools to get some metrics on data & alignments
# cat "$BAM_FILE" |
# while read i
# do
    # file=$(basename "$i")

    # Fetch filename from the array
    file=$(ls $ALIGNEDFOLDER/*.sorted.bam | sed "${SLURM_ARRAY_TASK_ID}q;d" | xargs -n 1 basename)

    echo \n">>> Computing alignment metrics for $file <<<"\n
    java -jar $PICARD $ALIGN \
        R=$GENOMEFOLDER/$GENOME \
        I=$ALIGNEDFOLDER/$file \
        O=$METRICSFOLDER/${file}_alignment_metrics.txt

    echo \n">>> Computing insert size metrics for $file <<<"\n
    java -jar $PICARD $INSERT \
        I=$ALIGNEDFOLDER/$file \
        OUTPUT=$METRICSFOLDER/${file}_insert_size_metrics.txt \
        HISTOGRAM_FILE=$METRICSFOLDER/${file}_insert_size_histogram.pdf

    echo \n">>> Computing coverage metrics for $file <<<"\n
    java -jar $PICARD $COVERAGE \
        R=$GENOMEFOLDER/$GENOME \
        I=$ALIGNEDFOLDER/$file \
        OUTPUT=$METRICSFOLDER/${file}_collect_wgs_metrics.txt\
        CHART=$METRICSFOLDER/${file}_collect_wgs_metrics.pdf

    #echo "Computing coverage for $file"
    samtools depth -a $ALIGNEDFOLDER/$file | gzip - > $METRICSFOLDER/${file}_coverage.gz

    echo \n">>> DONE! <<<"\n
#done 2> "$LOG_FOLDER"/03_metrics_"$TIMESTAMP".log
