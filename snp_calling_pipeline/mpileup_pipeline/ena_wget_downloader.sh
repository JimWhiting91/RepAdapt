#!/bin/bash
#SBATCH -J ena_wget_download
#SBATCH -o %x.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=0-12:00:00
#SBATCH --account=def-yeaman
#SBATCH -D .
##SBATCH --mail-type=ALL
##SBATCH --mail-user=james.whiting@ucalgary.ca

### Script takes RepAdapt metadata and downloads all fastq files and checks md5
RUN_ACC=$1
METADATA=$2
RAW_DATA=04_raw_data

# Download fastq to raw to raw_data dir
mkdir $RAW_DATA

# Fetch ftp addresses and md5
R1_FTP=$(grep $RUN_ACC $METADATA | cut -f15)
R2_FTP=$(grep $RUN_ACC $METADATA | cut -f16)
R1_MD5=$(grep $RUN_ACC $METADATA | cut -f17)
R2_MD5=$(grep $RUN_ACC $METADATA | cut -f18)

# Download each
wget $R1_FTP
wget $R2_FTP

# Move to raw data dir
mv $(basename $R1_FTP) $RAW_DATA/${RUN_ACC}_R1.fastq.gz
mv $(basename $R2_FTP) $RAW_DATA/${RUN_ACC}_R2.fastq.gz

# Check md5
md5_check1=($(md5sum $RAW_DATA/${RUN_ACC}_R1.fastq.gz))
md5_check2=($(md5sum $RAW_DATA/${RUN_ACC}_R2.fastq.gz))
if [ $md5_check1 == $R1_MD5 ]
then
  echo "$RUN_ACC R1 OK" >> $RAW_DATA/md5_checks.txt
else
  echo "$RUN_ACC R1 FAIL" >> $RAW_DATA/md5_checks.txt
fi

if [ $md5_check2 == $R2_MD5 ]
then
  echo "$RUN_ACC R2 OK" >> $RAW_DATA/md5_checks.txt
else
  echo "$RUN_ACC R2 FAIL" >> $RAW_DATA/md5_checks.txt
fi
