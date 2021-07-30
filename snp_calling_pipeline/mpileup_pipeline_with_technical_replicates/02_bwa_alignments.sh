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

# If this is our first run, make a list of all the trimmed reads for cleaning
if [ $SLURM_ARRAY_TASK_ID -eq 1 ]
then
  ls -1 $RAWDATAFOLDER/*R1.trimmed.fastq.gz | xargs -n 1 basename | sed 's/.R1.trimmed.fastq.gz//g' > $RAWDATAFOLDER/all_trimmed_ids.txt
fi

# Pull individual ID from the batch array
name=$(cut -f1 02_info_files/datatable.txt | sed "${SLURM_ARRAY_TASK_ID}q;d")

    # Name of uncompressed file
    file1=${name}.R1.trimmed.fastq.gz
    file2=${name}.R2.trimmed.fastq.gz
    echo ">>> Aligning file $file1 $file2 <<<
        "

    # Now clean if we have to
    # First check whether we need to edit the header
  if [ $(zcat $RAWDATAFOLDER/$file1 | head -n100 | grep "@" | cut -d" " -f1 | sed 's/\./\t/g' |  awk '{print $NF}' | sort | uniq | wc -l) -eq 1 ];
  then

    echo ">>> R1 and R2 files for $name have malformed headers. Cleaning and removing suffixes.
    "

    mv $RAWDATAFOLDER/$file1 $RAWDATAFOLDER/${name}.R1.trimmed_dirty.fastq.gz
    mv $RAWDATAFOLDER/$file2 $RAWDATAFOLDER/${name}.R2.trimmed_dirty.fastq.gz
    zcat $RAWDATAFOLDER/${name}.R1.trimmed_dirty.fastq.gz | sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" | gzip > $RAWDATAFOLDER/$file1
    zcat $RAWDATAFOLDER/${name}.R2.trimmed_dirty.fastq.gz | sed -E "s/^((@|\+)SRR[^.]+\.[^.]+)\.(1|2)/\1/" | gzip > $RAWDATAFOLDER/$file2

    # Keep clean
    rm -f $RAWDATAFOLDER/${name}.R1.trimmed_dirty.fastq.gz $RAWDATAFOLDER/${name}.R2.trimmed_dirty.fastq.gz

  else
    echo ">>> R1 and R2 files for $name passed QC
    "
  fi

    #statements

    # Set ID
    ID="@RG\tID:ind\tSM:ind\tPL:Illumina"

    # echo "$name"
    # echo "$name2"

    # Align reads
    bwa mem -t $NCPU -R $ID $GENOMEFOLDER/$GENOME $RAWDATAFOLDER/$file1 $RAWDATAFOLDER/$file2 |
    samtools view -Sb -q 10 - > $ALIGNEDFOLDER/${name%}.bam

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
