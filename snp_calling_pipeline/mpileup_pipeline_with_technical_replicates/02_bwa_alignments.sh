#!/bin/bash

#SBATCH -J 02.BWA
#SBATCH -o 98_log_files/%x_%A_array%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=0-12:00:00

# Load needed modules
module load bwa samtools

#cd $SLURM_SUBMIT_DIR

##Keep some info. about the run/script
TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)
SCRIPT=$0
NAME=$(basename $0)
LOG_FOLDER="98_log_files"
cp $SCRIPT $LOG_FOLDER/"$TIMESTAMP"_"$NAME"_%J

# Global variables
GENOMEFOLDER="03_genome"
GENOME=$(ls -1 $GENOMEFOLDER/*{fasta,fa,fasta.gz,fa.gz} | xargs -n 1 basename)
INDGENOME=${GENOME}.fai
RAWDATAFOLDER="05_trimmed_data"
ALIGNEDFOLDER="06_bam_files"
#SAMPLE_FILE="02_info_files/samples_split/bwa0000"

TIMESTAMP=$(date +%Y-%m-%d_%Hh%Mm%Ss)

# Test if user specified a number of CPUs
if [[ -z "$NCPU" ]]
then
    NCPU=4
fi

# # Iterate over sequence file pairs and map with bwa
# cat "$SAMPLE_FILE" |
# while read file
# do

# Pull individual ID from the batch array
name=$(ls -1 $RAWDATAFOLDER/*1.trimmed.fastq.gz | xargs -n 1 basename | sed 's/.R1.trimmed.fastq.gz//g' | sed "${SLURM_ARRAY_TASK_ID}q;d")

    # Name of uncompressed file
    file1=${name}.R1.trimmed.fastq.gz
    file2=${name}.R2.trimmed.fastq.gz
    echo ">>> Aligning file $file1 $file2 <<<
        "

    # #######################################
    # # Uncomment these lines in rare cases where SRA download has suffixed R1 and R2 fastq headers with .1 and .2, which errors out bwa. This removes them
    # ###
    # mv $RAWDATAFOLDER/$file1 $RAWDATAFOLDER/${name}.R1.trimmed_dirty.fastq.gz
    # mv $RAWDATAFOLDER/$file2 $RAWDATAFOLDER/${name}.R2.trimmed_dirty.fastq.gz
    # zcat $RAWDATAFOLDER/${name}.R1.trimmed_dirty.fastq.gz | sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" | gzip > $RAWDATAFOLDER/$file1
    # zcat $RAWDATAFOLDER/${name}.R2.trimmed_dirty.fastq.gz | sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" | gzip > $RAWDATAFOLDER/$file2
    # #######################################

    # Set ID
    ID="@RG\tID:ind\tSM:ind\tPL:Illumina"

    # echo "$name"
    # echo "$name2"

    # Align reads
    bwa mem -t $NCPU -R $ID $GENOMEFOLDER/$GENOME $RAWDATAFOLDER/$file1 $RAWDATAFOLDER/$file2 |
    samtools view -Sb -q 10 - > $ALIGNEDFOLDER/${name%}.bam

    bwa mem -t $NCPU -R $ID $GENOMEFOLDER/$GENOME $RAWDATAFOLDER/test1.fastq.gz $RAWDATAFOLDER/test1.fastq.gz |
    samtools view -Sb -q 10 - > $ALIGNEDFOLDER/test.bam

    # Sort
    samtools sort --threads $NCPU $ALIGNEDFOLDER/${name%.R1.trimmed.fastq.gz}.bam \
        > $ALIGNEDFOLDER/${name%.R1.trimmed.fastq.gz}.sorted.bam

    # Index
    samtools index $ALIGNEDFOLDER/${name%.R1.trimmed.fastq.gz}.sorted.bam

    &> $LOG_FOLDER/02_mapping_${name}_${TIMESTAMP}.log
#done

echo " >>> Cleaning a bit...
"
#rm "$ALIGNEDFOLDER"/"${name%}".bam
echo "Completed SLURM job $SLURM_JOB_ID in $(sacct -nXj $SLURM_JOB_ID -o elapsed)"
