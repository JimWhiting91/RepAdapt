#!/bin/bash

# Set up directory and move to top of tree
MAIN=/lu213/james.whiting/RepAdapt
cd $MAIN

# Currently input is taken as the VCF
# Each VCF however may include multiple species, which are separated during analysis
#VCF=$1

# For now specify manually based on the full path
# VCF=/lu213/james.whiting/RepAdapt/data/VCFs/10_Capsella_Weigel/Co_full_concatened.vcf.gz
METADATA=metadata/sample_species_vcf_author_map_v2_210519.csv

# Set general variables
POP_STRUCTURE_SNPN=10000
NCORE=16
TC_THRESH=0.99

# Try and loop over these VCF
for VCF in $(grep -v "ool" $METADATA | sed 's/,/\t/g' | cut -f10 | uniq | grep "vcf" | grep -v "Murray" | grep -v "Capsella")
do

  echo "STARTING TEST GEA FOR $VCF"

  # Fetch the reference genome from the VCF header, and also the GFF
  REFERENCE=$(gunzip -c $VCF | grep -m1 "##reference")
  REFERENCE=$(echo ${REFERENCE##*/})
  REFERENCE_PATH=$(find data/reference_genomes/* -name $REFERENCE)
  REFERENCE_DIR=$(echo "${REFERENCE_PATH%/*}")

  # This currently struggles with gzipped gffs...
  GFF=$(ls -d $REFERENCE_DIR/*.gff* | grep -v "gz")

  # Write these to log file
  echo "LOGFILE FOR GEA OF $VCF" > $VCF.gea.log
  echo "REFERENCE = $REFERENCE" >> $VCF.gea.log
  echo "GFF = $GFF" >> $VCF.gea.log

  # First, quantify population structure
  bin/Rscript R/quantify_population_structure.R $VCF $METADATA $NCORE $POP_STRUCTURE_SNPN >> $VCF.gea.log

  # Run the GEAs over climate data
  # Something in here to define which climate data to use?
  bin/Rscript R/perform_GEA.R $VCF $METADATA $NCORE >> $VCF.gea.log
  # Or, using argparse if possible
  #bin/Rscript perform_GEA.R -v $VCF -m $METADATA -n 24

  # Collate results to gene-based WZA
  bin/Rscript R/collate_GEA_to_gene_WZA.R $VCF $METADATA $NCORE $GFF $TC_THRESH >> $VCF.gea.log

  # Add something here to tarball up .frq files as they are bloody big.

done
