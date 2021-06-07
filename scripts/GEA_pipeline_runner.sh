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

# Loop over these
pool_test=(/lu213/james.whiting/RepAdapt/data/VCFs/02_Alyrata_Willi_pool/Alyr_full_concatened.vcf.gz /lu213/james.whiting/RepAdapt/data/VCFs/11_Athaliana_Gunther_pool/Athal_full_concatened.vcf.gz /lu213/james.whiting/RepAdapt/data/VCFs/12_Athaliana_Roux_pool/Athal_full_concatened.vcf.gz)
pool_test_ref=(Alyr_v.1.0_genomic.fasta GCF_000001735.4_TAIR10.1_genomic.fasta GCF_000001735.4_TAIR10.1_genomic.fasta)

for i in {0..2}
do
# one off
VCF="${pool_test[$i]}"
REFERENCE="${pool_test_ref[$i]}"
IS_POOL=true

# Try and loop over these VCF
#for VCF in $(grep -v "ool" $METADATA | sed 's/,/\t/g' | cut -f10 | uniq | grep "vcf" | grep -v "Murray" | grep -v "Capsella")
# for VCF in $(grep -v "ool" $METADATA | sed 's/,/\t/g' | cut -f10 | uniq | grep "vcf" | grep -v "Murray" | grep -v "Capsella" | grep -v "sylvestris" | grep -v "bursa" | grep -v "tremula")
# do

  echo "STARTING TEST GEA FOR $VCF"

  # Fetch the reference genome from the VCF header, and also the GFF
  # REFERENCE=$(gunzip -c $VCF | grep -m1 "##reference")
  # REFERENCE=$(echo ${REFERENCE##*/})
  REFERENCE_PATH=$(find data/reference_genomes/* -name $REFERENCE)
  REFERENCE_DIR=$(echo "${REFERENCE_PATH%/*}")

  # This currently struggles with gzipped gffs...
  GFF=$(ls -d $REFERENCE_DIR/*.gff* | grep -v "gz")

  # Write these to log file
  echo "LOGFILE FOR GEA OF $VCF" > $VCF.gea.log
  echo "REFERENCE = $REFERENCE" >> $VCF.gea.log
  echo "GFF = $GFF" >> $VCF.gea.log

  ################################################################################################################
  # STEP 1: Summarise Population Structure and Dataset Features

#   # First, quantify population structure
# if [ "IS_POOL"=true ]
# then
#   bin/Rscript R/quantify_population_structure.R \
#   --vcf=$VCF \
#   --metadata=$METADATA \
#   --n_cores=$NCORE \
#   --sub_SNP=$POP_STRUCTURE_SNPN \
#   --pool >> $VCF.gea.log
# else
#   bin/Rscript R/quantify_population_structure.R \
#   --vcf=$VCF \
#   --metadata=$METADATA \
#   --n_cores=$NCORE \
#   --sub_SNP=$POP_STRUCTURE_SNPN >> $VCF.gea.log
# fi

  ################################################################################################################
  # STEP 2: GEA over worldclim

  # Run the GEAs over climate data - Make sure to set IS_POOL to organise how AFs are handled
  # Something in here to define which climate data to use?
if [ "IS_POOL"=true ]
then
  bin/Rscript R/perform_GEA.R \
    --vcf=$VCF \
    --metadata=$METADATA \
    --n_cores=$NCORE \
    --pool >> $VCF.gea.log
else
  bin/Rscript R/perform_GEA.R \
    --vcf=$VCF \
    --metadata=$METADATA \
    --n_cores=$NCORE >> $VCF.gea.log
fi

  ################################################################################################################
  # STEP 3: Summarise GEA to WZA and TC

  # Collate results to gene-based WZA
  bin/Rscript R/collate_GEA_to_gene_WZA.R $VCF $METADATA $NCORE $GFF $TC_THRESH >> $VCF.gea.log

  # Add something here to tarball up .frq files as they are bloody big.

done
