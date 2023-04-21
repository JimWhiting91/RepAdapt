#!/bin/bash

#SBATCH -J "05b.merge"
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=50G
#SBATCH --time=00-01:30:00

#cd $SLURM_SUBMIT_DIR

#----- This script performs merging across bams that represent mulitple sequence files from one individual
#----- Scrpt is based on script 5 here: https://github.com/josieparis/gatk-snp-calling

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"
cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$NAME"


# Load needed modules - ComputeCanada clusters
module load picard java

export JAVA_TOOL_OPTIONS="-Xms2g -Xmx50g "
export _JAVA_OPTIONS="-Xms2g -Xmx50g "

# Global variables
INBAM="06_bam_files"
OUTBAM="06_bam_files"
ADDRG="AddOrReplaceReadGroups"
PICARD=$EBROOTPICARD/picard.jar
DATATABLE=02_info_files/datatable.txt

# Remove duplicates from bam alignments
echo "Merging samples...
"

# Fetch filename from the array
SAMPLE=$(cut -f14 $DATATABLE | sort | uniq | sed "${SLURM_ARRAY_TASK_ID}q;d")
READ_ID_ARRAY=( $(awk -F "\t" -v BAM=$SAMPLE '$14==BAM {print $1}' $DATATABLE) )

# Find all bams for this individual and merge, checking first that the RG file does actually exist
cmd=""
for READ_ID in "${READ_ID_ARRAY[@]}"
do
	rg_bam=${OUTBAM}/${READ_ID}_RG.bam
	if test -f "$rg_bam"; then
	cmd+="I=${OUTBAM}/${READ_ID}_RG.bam "
fi
done

# Now merge
java -Xmx10g -jar $PICARD MergeSamFiles $cmd O=$OUTBAM/${SAMPLE}.merged.bam TMP_DIR=~/scratch/

echo " >>> Cleaning a bit...
"
#rm "$INBAM"/"$file"
echo "
DONE! Check your files"