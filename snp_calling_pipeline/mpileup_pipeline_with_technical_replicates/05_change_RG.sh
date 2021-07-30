#!/bin/bash

#SBATCH -J "05.RG"
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=7G
#SBATCH --time=00-01:30:00

#cd $SLURM_SUBMIT_DIR

# Copy script to log folder
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"
cp "$SCRIPT" "$LOG_FOLDER"/"$TIMESTAMP"_"$NAME"

# Load needed modules - ComputeCanada clusters
#module load nixpkgs/16.09
module load picard java samtools
#module load picard/2.20.6


# Global variables
INBAM="06_bam_files"
OUTBAM="06_bam_files"
ADDRG="AddOrReplaceReadGroups"
PICARD=$EBROOTPICARD/picard.jar
DATATABLE=02_info_files/datatable.txt
#DEDUP_FILE="02_info_files/dedup_split/dedup0000"

# Remove duplicates from bam alignments
echo "Editing RG...
"
# cat "$DEDUP_FILE" |
# while read file
# do

# Fetch filename from the array
sample_name=$(cut -f1 02_info_files/datatable.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")
file=${sample_name}.dedup.bam

# Fetch all our RG info...
#new_RGID=$(grep $sample_name $DATATABLE | cut -f10)
new_RGSM=$(grep $sample_name $DATATABLE | cut -f14)

        echo "
             >>> Computing RG for $file<<<
             "
        java -jar $PICARD $ADDRG \
	    I=$INBAM/$file \
	    O=$OUTBAM/${sample_name}_RG.bam \
	    RGID=${sample_name} \
	    RGLB=${sample_name}_LB \
	    RGPL=ILLUMINA \
	    RGPU=unit1 \
	    RGSM=${new_RGSM}
        # Index
        echo "
            >>> Indexing ${sample_name}_RG.bam <<<
            "
        samtools index $INBAM/${sample_name}_RG.bam
# done 2> "$LOG_FOLDER"/05_RG_"$TIMESTAMP".log

echo " >>> Cleaning a bit...
"
#rm "$INBAM"/"$file"
echo "
DONE! Check your files"

#echo "
#I also cleaned out intermediate files. Cheers!"
#rm $INBAM/*.sorted.bam && rm $INBAM/*.sorted.bam.bai
