#!/bin/bash

# Script uses agat gff tools to build a proteome for a given reference genome with annotation
GFF=$1
GENOME=$2
OUT_NAME="${GFF%.*}"

# Load the conda environment
#conda activate agat-env

# Extract the longest isoform and save this to a temporary new gff
agat_sp_keep_longest_isoform.pl --gff $GFF -o ${OUT_NAME}_long_iso.gff

# Fetch the protein sequences
agat_sp_extract_sequences.pl --gff ${OUT_NAME}_long_iso.gff -f $GENOME -p -o ${OUT_NAME}_AA_proteome.fa

# Clean
rm -f ${OUT_NAME}_long_iso.gff 
