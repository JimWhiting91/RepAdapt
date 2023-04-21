#!/bin/bash
#SBATCH -J sra_download
#SBATCH -o %x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=0-12:00:00
#SBATCH --account=def-yeaman
#SBATCH -D .

# module load nixpkgs/16.09 sra-toolkit/2.9.6
module load StdEnv/2020 gcc/9.3.0 sra-toolkit/3.0.0

# SRR... code should given as a trailing variable...
SRR_CODE=$1

# Download sra to raw
mkdir 04_raw_data
cd 04_raw_data/
# And get fastqs
mkdir fastq_files

# Check whether already exists
prefetch -X 16G -T fastq -O ./ $SRR_CODE

# Need to catch an exception for cases when the download is downloaded via ENA
if [[ $(echo $SRR_CODE | head -c 3) -eq "ERR" ]]
then
  # Change the names to what the metadata is expecting
  mv $SRR_CODE/*R1*gz fastq_files/${SRR_CODE}_R1.fastq.gz
  mv $SRR_CODE/*R2*gz fastq_files/${SRR_CODE}_R2.fastq.gz

else
  # Validate these
  echo ">>> Validating $SRR_CODE.sra
  "
  vdb-validate $SRR_CODE.sra

  # Dump fastqs
  #for file in $(ls *.sra)
  #do
  fasterq-dump -S -v -O ./fastq_files/ --threads 4 $SRR_CODE.sra
  #done

  # Gzip them all
  cd fastq_files
  parallel -j4 "gzip {}" ::: $(ls ${SRR_CODE}*)

  # Rename stupid sra naming
  for fastq in $(ls ${SRR_CODE}.*.gz)
  do
   mv $fastq $(echo $fastq | sed 's/.sra_/_R/')
  done

fi

cd ..
