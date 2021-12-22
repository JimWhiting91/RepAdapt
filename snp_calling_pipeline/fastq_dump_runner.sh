#!/bin/bash
#SBATCH -J fastq_dump_runner_kremer
#SBATCH -o %x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=0-168:00:00
#SBATCH --account=def-yeaman
#SBATCH --mail-type=ALL
#SBATCH --mail-user=james.whiting@ucalgary.ca

cd /home/jimw91/scratch/pool_snpcalling/kremer_Qpetraea

module load nixpkgs/16.09 sra-toolkit/2.9.6

# Fetch the .sra files
SRA_FILES=""
for file in $(cat 02_info_files/sra_download_list.txt)
do
SRA_FILES+="$file "
done

# Download sra to raw
mkdir 04_raw_data
cd 04_raw_data/
prefetch -X 100G -T fastq -O ./ $SRA_FILES

# Validate these
for file in $(ls *.sra)
do
echo ">>> Validating $file
"
vdb-validate $file
done

# And get fastqs
mkdir fastq_files

# Dump fastqs
for file in $(ls *.sra)
do
fasterq-dump -S -v -O ./fastq_files/ --threads 6 $file
done

# Gzip them all
cd fastq_files
parallel -j6 "gzip {}" ::: $(ls *fastq)
