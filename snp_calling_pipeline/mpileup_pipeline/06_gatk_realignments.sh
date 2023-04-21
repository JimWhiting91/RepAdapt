#!/bin/bash
# 1 CPU
# 30 Go

#SBATCH -J "06.ReAlignmentS"
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=00-168:00:00

#cd $SLURM_SUBMIT_DIR

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"
cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$NAME"

# Load needed modules
module load StdEnv/2020 samtools/1.12

# Uncomment these for big genomes
export JAVA_TOOL_OPTIONS="-Xms2g -Xmx32g "
export _JAVA_OPTIONS="-Xms2g -Xmx32g "

# Global variables
BAM="06_bam_files"
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/*{fasta,fa,fasta.gz,fa.gz} | xargs -n 1 basename)
INDGENOME=${GENOME}.fai
DATATABLE=02_info_files/datatable.txt

# Build Bam Index
echo " >>> Realigning...
"


# Fetch filename from the array
sample_name=$(cut -f14 $DATATABLE | sort | uniq | sed "${SLURM_ARRAY_TASK_ID}q;d")
file=${sample_name}.merged.bam

    echo "
         >>> Realigning TARGET for $file <<<
         "

    # Index the bam file First
    samtools index $BAM/$file

    # Now load modules
    module purge
    module load nixpkgs/16.09
    module load java gatk/3.8

    # Realign
    java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
        -T RealignerTargetCreator \
        -R $GENOMEFOLDER/$GENOME \
        -I $BAM/$file \
        -o $BAM/${sample_name}.intervals

    echo "
         >>> Realigning INDELs for $file <<<
         "
    java -jar $EBROOTGATK/GenomeAnalysisTK.jar \
        -T IndelRealigner \
        -R $GENOMEFOLDER/$GENOME \
        -I $BAM/$file \
        -targetIntervals $BAM/${sample_name}.intervals \
        --consensusDeterminationModel USE_READS  \
        -o $BAM/${sample_name}.realigned.bam

echo ">>> Cleaning a bit...
"
echo "
DONE! Check your files"