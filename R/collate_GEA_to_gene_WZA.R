####################################################################################
# This script takes a VCF input and performs GEA analysis over climate variables
lib <- c("parallel","data.table","sp","VGAM","Rfast")
sapply(lib,library,character.only=T)

# Read the VCF in from the command line
args <- commandArgs(TRUE)
vcf_path <- as.character(args[1])
gff_path <- as.character(args[2])
metadata_path <- as.character(args[2])
n_cores <- as.character(args[3])

# For setting up purposes
vcf_path <- "/lu213/james.whiting/RepAdapt/data/VCFs/06_Ahalleri_Kubota_wgs/Ahalg_full_concatened.vcf.gz"
results_dir <- "outputs/Arabidopsis_halleri_Kubota_Individual/"
gff_path <- "data/reference_genomes/A.halleri_gemmifera_reference/Arabidopsis_halleri.Ahal2.2.45.gff3"
metadata_path <- "metadata/sample_species_vcf_author_map_v2_210507.csv"

n_cores <- 16

################################################
# Function library

################################################
##### SET UP CLIMATE VARS #####

climate_vars <- c("mean_temp","mean_diurnal","isothermality","temp_seasonality","max_temp_warmest_month",
                        "min_temp_coldest_month","temp_range","mean_temp_wet_quarter","mean_temp_dry_quarter","mean_temp_warm_quarter",
                        "mean_temp_cold_quarter","annual_precip","precip_wet_month","precip_dry_month","precip_seasonality","precip_wet_quarter",
                        "precip_dry_quarter","precip_warm_quarter","precip_cold_quarter")
names(climate_vars) <- c("Annual Mean Temperature",
                  "Mean Diurnal Range",
                  "Isothermality",
                  "Temperature Seasonality",
                  "Max Temperature of Warmest Month",
                  "Min Temperature of Coldest Month",
                  "Temperature Annual Range",
                  "Mean Temperature of Wettest Quarter",
                  "Mean Temperature of Driest Quarter",
                  "Mean Temperature of Warmest Quarter",
                  "Mean Temperature of Coldest Quarter",
                  "Annual Precipitation",
                  "Precipitation of Wettest Month",
                  "Precipitation of Driest Month",
                  "Precipitation Seasonality",
                  "Precipitation of Wettest Quarter",
                  "Precipitation of Driest Quarter",
                  "Precipitation of Warmest Quarter",
                  "Precipitation of Coldest Quarter")

################################################
##### PREPARE GENE SNP OVERLAPS #####
# Set gene boundaries from gff
gff <- read.table(gff_path,fill=TRUE,comment.char="#",sep="\t")
colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","att")

# Set up bed file for all genes
gene_dd <- data.frame(chr=gff[gff$type=="gene","seqid"],
                        start=gff[gff$type=="gene","start"],
                        end=gff[gff$type=="gene","end"],
                        info=gff[gff$type=="gene","att"])

# Set up keyed index for overlap
gene_bed <- data.table(gene_dd[,1:3])
setkey(gene_bed, chr, start, end)


# Set up SNP boundaries
GEA_res <- list.files(results_dir,pattern="GEA.tsv")
GEA_SNPs <- data.frame(fread(paste0(results_dir,"/",GEA_res[1]),fill=T))[,1]

# Fetch chromosomes and SNPs
snp_dd <- data.frame(chr=sapply(strsplit(GEA_SNPs,":"),'[[',1),
                     start=as.integer(sapply(strsplit(GEA_SNPs,":"),'[[',2)))

# If a species identifier has been added, remove it
if(any(grep("_",snp_dd[,1]))){
  snp_dd$chr <- sapply(strsplit(snp_dd$chr,"_"),'[[',2)
}

# Set these up as genomic regions
snp_dd$end <- snp_dd$start+1
snp_bed <- data.table(snp_dd)
setkey(snp_bed, chr, start, end)

# Find overlaps
gene_snp_overlap <- data.frame(foverlaps(gene_bed, snp_bed, minoverlap=1, nomatch = NA))

# Use overlaps to group genes
gene_snp_overlap$snp_id <- paste0(gene_snp_overlap$chr,":",gene_snp_overlap$start)
gene_snp_overlap$gene_id <- paste0(gene_snp_overlap$chr,":",gene_snp_overlap$i.start,"-",gene_snp_overlap$i.end)

# And add these to original beds
gene_dd$gene_id <- paste0(gene_dd$chr,":",gene_dd$start,"-",gene_dd$end)
snp_dd$snp_id <- paste0(snp_dd$chr,":",snp_dd$start)

# Make a list of SNP IDs for relevant genes...
gene_snp_list <- pbmcapply::pbmclapply(unique(gene_snp_overlap$gene_id), function(x) return(gene_snp_overlap[gene_snp_overlap$gene_id == x,"snp_id"]),mc.cores=n_cores)

################################################
##### PREPARE ALLELE FREQUENCY INPUT #####

# Based on the VCF, read in the metadata and subset it
metadata <- read.csv(metadata_path)
sub_meta <- metadata[metadata$VCF==vcf_path,]

# Get the individuals we're getting frequencies for
inds <- sub_meta[,"Sample_name"]
inds_to_keep <- paste0("--indv ",inds,collapse = " ")

# Calculate global allele frequencies with vcftools only keeping individuals with metadata
system(paste0("vcftools --gzvcf ",vcf_path," --freq ",inds_to_keep," --out ",results_dir,"/whole_data_AFs"),wait = T)

# Fetch the allele frequency results
freq_res <- data.frame(fread(paste0(results_dir,"/whole_data_AFs.frq"),fill=T))
colnames(freq_res) <- c("chr","bp","allele_N","allele_count","ref_freq","alt_freq")

# Format for analysis
freq_dd <- data.frame(chr=freq_res$chr,
                       bp=freq_res$bp,
                       ref_freq=as.numeric(sapply(strsplit(freq_res$ref_freq,":"),'[[',2)),
                       alt_freq=as.numeric(sapply(strsplit(freq_res$alt_freq,":"),'[[',2)))

# And make a vector of minor allele frequencies
freq_dd$minor_af <- rje::rowMins(freq_dd[,c("ref_freq","alt_freq")])

################################################
##### ESTIMATE WZA FOR EACH GENE #####
## This code assumes that you have a dataframe (called GEA)
## Each row corresponds to a single SNP
##
## The dataframe should have at least three columns:
##
## GEA$pVal - the p-value for a particular SNP - could be spearman's rho or whatever you happen to habe
## GEA$pbar_qbar - the average minor allele frequency across populations (pbar) multiplied by the average major allele frequency
## GEA$gene - the name of the gene or genomic region that a SNP corresponds to

## Convert one-sided p-values to Z-scores
GEA$z <- qnorm(GEA$pVal , lower.tail = F)
## Calculate the numerator of the Weighted-Z score
weiZ_num <- tapply( GEA$pbar_qbar * GEA$z, GEA$gene, sum )
## Calculate the denominator of the Weighted-Z score
weiZ_den <- sqrt(tapply( GEA$pbar_qbar^2, GEA$gene, sum ))
## Bring data together into a new DF
Z_df <- data.frame( gene = names(weiZ_num), weiZ = weiZ_num/weiZ_den)
