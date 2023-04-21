#!/bin/bash

# Set our working directory
MAIN=~/RepAdapt/
cd $MAIN

# Set map file for genomes and their annotations
GENOME_GFF_MAP=metadata/vcf_genome_gff_210805_map.txt

# Where are our genomes?
GENOME_DIR=data/reference_genomes
PROTEOME_DIR=data/proteomes
RUN_NAME=221213_18_genomes_Ptaeda_isoforms_removed
NUMCORE=16
AA_FILTER=50

# Set up our conda env
conda activate orthology-env

# Make all proteomes if needed
grep -vw "NA" $GENOME_GFF_MAP | grep -v "#" | cut -f2,3 | sort -nk1 | uniq > tmp.map
for GFF in $(cut -f2 tmp.map)
do
  GENOME=$(grep $GFF tmp.map | cut -f1)
  echo "PREPARING PROTEOME FOR $GENOME
  "
  bash scripts/processing/prepare_proteome.sh $GFF $GENOME
done

# Fetch all the proteomes to a file and copy them over to somewhere new
find $GENOME_DIR/* -name "*proteome.fa" > metadata/proteome_paths.txt
mkdir $PROTEOME_DIR $PROTEOME_DIR/filtered
for prot in $(cat metadata/proteome_paths.txt)
do
  cp $prot data/proteomes/
done

# Rename these for convenience
for prot in $(ls data/proteomes/)
do
  GFF_NAME=$(echo $prot | sed 's/_AA_proteome.fa//g')
  SPECIES_NAME=$(grep $GFF_NAME $GENOME_GFF_MAP | cut -f5 | head -n1)
  mv data/proteomes/$prot data/proteomes/${SPECIES_NAME}.faa

  # Also remove proteins with less than AA_FILTER residues
  bin/seqkit seq -m $AA_FILTER $PROTEOME_DIR/${SPECIES_NAME}.faa > $PROTEOME_DIR/filtered/${SPECIES_NAME}.faa
done

# Run OrthoFinder
orthofinder -f $PROTEOME_DIR \
-o outputs/orthology/orthofinder_res_${RUN_NAME} \
-n $RUN_NAME \
-t $NUMCORE \
-a $NUMCORE

conda deactivate
