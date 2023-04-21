#!/bin/bash

# Set up directory and move to top of tree
MAIN=/lu213/james.whiting/RepAdapt
cd $MAIN
PIPELINE_DIR=R/01_GEA_pipeline

# Set a run_name...
RUN_NAME=230321

# For now specify manually based on the full path
VCF_GENOME_GFF=metadata/vcf_genome_gff_220830_map.txt

# Make a file of VCFs to run if need be
# This should just be a simple text file listing relative paths from $MAIN to VCFs. VCFs can be commented out if already run.
VCF_TO_RUN=metadata/vcfs_to_run.txt

# Set general variables
POP_STRUCTURE_SNPN=10000  # How many SNPs to use for popstructure analyses
WZA_DOWNSAMPLE=0.75 # How many SNPs per gene to keep for WZA downsampling. If <1 we take it as a quantile of the distribution of per gene snp density...
GENE_FLANK_SIZE=500 # What size flanking regions to add onto genes when collating SNPs to gene scores...
NCORE=12
TC_THRESH=0.99 # What upper % for dbinom with top-candidate approach
RDA_DOWNSAMPLE=10000  # How many SNPs to downsample for RDA analyses

# Loop over VCF that we are ready to run and that haven't been run
for VCF in $(grep -v "#" $VCF_TO_RUN)
do

  # Set Dataset Variables...
  REFERENCE=$(grep $VCF $VCF_GENOME_GFF | awk '{print $2}')
  GFF=$(grep $VCF $VCF_GENOME_GFF | awk '{print $3}')
  IS_POOL=$(grep $VCF $VCF_GENOME_GFF | awk '{print $4}')
  SPECIES_CODE=${RUN_NAME}_$(grep $VCF $VCF_GENOME_GFF | awk '{print $8}')

  # Write these to a log file
  echo "Run name = $RUN_NAME" > $VCF.gea.log
  echo "ref_genome = ${REFERENCE}" >> $VCF.gea.log
  echo "gff = ${GFF}" >> $VCF.gea.log
  echo "is_pool = ${IS_POOL}" >> $VCF.gea.log

  ################################################################################################################
  # STEP 1: GEA over worldclim
  echo ">>> STARTING ALL GEA FOR $VCF
  "

  # Run the GEAs over climate data - Make sure to set IS_POOL to organise how AFs are handle
  if [ "$IS_POOL" = "TRUE" ]
  then
    bin/Rscript $PIPELINE_DIR/perform_GEA.R \
      --dataset_dir=$SPECIES_CODE \
      --vcf=$VCF \
      --n_cores=$NCORE \
      --pool >> $VCF.gea.log
  else
    bin/Rscript $PIPELINE_DIR/perform_GEA.R \
      --dataset_dir=$SPECIES_CODE \
      --vcf=$VCF \
      --n_cores=$NCORE >> $VCF.gea.log
  fi

  ################################################################################################################
  # STEP 2: Summarise Population Structure and Dataset Features
  echo ">>> STARTING POP STRUCTURE SUMMARIES FOR $VCF
  "
  #  Secondly, quantify population structure
  if [ "$IS_POOL" = "TRUE" ]
  then
    bin/Rscript $PIPELINE_DIR/quantify_population_structure.R \
    --dataset_dir=$SPECIES_CODE \
    --vcf=$VCF \
    --n_cores=$NCORE \
    --sub_SNP=$POP_STRUCTURE_SNPN \
    --pool >> $VCF.gea.log
  else
    bin/Rscript $PIPELINE_DIR/quantify_population_structure.R \
    --dataset_dir=$SPECIES_CODE \
    --vcf=$VCF \
    --n_cores=$NCORE \
    --sub_SNP=$POP_STRUCTURE_SNPN >> $VCF.gea.log
  fi

  ###############################################################################################################
  # STEP 3: Summarise GEA to WZA and TC
  echo ">>> STARTING GEA > WZA SUMMARIES FOR $VCF
  "
  # Collate results to gene-based WZA
    bin/Rscript $PIPELINE_DIR/collate_GEA_to_gene_WZA_empirical_pvals.R \
    --dataset_dir=$SPECIES_CODE \
    --vcf=$VCF \
    --n_cores=$NCORE \
    --gff=$GFF \
    --snp_per_gene=$WZA_DOWNSAMPLE \
    --gene_flank_size=$GENE_FLANK_SIZE \
    --tc_threshold=$TC_THRESH >> $VCF.gea.log

  ################################################################################################################
  # STEP 4: Partition genetic variance among climate variables
  echo ">>> STARTING RDA VARIANCE PARTITIONING FOR $VCF
  "
    bin/Rscript $PIPELINE_DIR/partition_variance_with_RDA.R \
    --vcf=$VCF \
    --dataset_dir=$SPECIES_CODE \
    --snp_downsample=$RDA_DOWNSAMPLE \
    --n_cores=$NCORE >> $VCF.gea.log

  # Finish loop
  echo ">>> Full GEA pipeline finished for $VCF
  ################################################################################################################
  ################################################################################################################
  "
done
