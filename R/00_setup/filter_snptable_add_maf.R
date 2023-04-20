################################################################################################################
# Filters SNP Table based on Brandon Lind recommendations and calculates global allele frequency for all sites
# Load these
lib <- c("parallel","data.table","sp","VGAM","argparse","R.utils","tidyverse")
dummy <- suppressPackageStartupMessages(lapply(lib,function(x){ 
  if (!require(x,character.only = T)){ 
    install.packages(x,repos='https://utstat.toronto.edu/cran/')
    library(x,character.only = T)
  } else {
    library(x,character.only = T)
  }
}))

# Filters as described:
# Mask sites where GQ < 20
# Filter for <25% missing data across pools
# MAF >= 0.05
# Remove sites from soft-masked reference
# INDELS are already removed by bcftools prior to filtering

# Filtering...
GQ_threshold = 20
max_missing_threshold = 0.25
maf_threshold = 0.05

# Set up environment ------------------------------------------------------
args <- commandArgs(asValues = T,excludeReserved = T)[-1]
print(args)

snptable_path <- args$snptable
metadata_path <- args$metadata

if(is.null(snptable_path) | is.null(metadata_path)){
  stop("ERROR: No SNPTable and/or metadata provided")
}
# if(is.null(n_cores)){
#   n_cores <- 1
# }
# 
# # Set up parallel cores
# cl <- makeCluster(n_cores, type="FORK")
# registerDoParallel(cl)
# 
# # Kill cluster
# on.exit(stopCluster(cl))

# Start filtering ---------------------------------------------------------
# For now...
# snptable_path="data/VCFs/11_Athaliana_Gunther_pool/test_table.txt"
# n_cores=16

# Read in the file
snptable <- data.frame(fread(snptable_path))

# Fetch the names of our individual pools
pool_names <- gsub(".GT","",grep(".GT",colnames(snptable),value = T))

# Read in ploidy and filter if needed -------------------------------------
# Make a ploidy file from the metadata
metadata <- read.csv(metadata_path)
ploidy_file <- metadata[grep(dirname(snptable_path),metadata$VCF),c("Sample_name","Pool_Size")]
colnames(ploidy_file) <- c("pool","ploidy")
ploidy_file <- na.omit(ploidy_file[match(pool_names, ploidy_file[,1]),])


# Remove individuals if we need to
if(nrow(ploidy_file) < length(pool_names)){
  message(paste0(">>> Removing ",length(pool_names) - nrow(ploidy_file)," pool with missing metadata...
                 "))
  to_remove <- pool_names[!(pool_names %in% ploidy_file$pool)]
  to_remove_cols <- unlist(sapply(to_remove,function(x) grep(x,colnames(snptable))))
  
  # And remove..
  snptable <- snptable[,!(colnames(snptable) %in% to_remove_cols)]
  pool_names <- pool_names[!(pool_names %in% to_remove)]
}

##### Run through each pool and mask any cases where GQ < 20 ##### 
cat(paste0(">>> Masking based on GQ of ",GQ_threshold,"\n"))
for(pool in pool_names){
  
  # For data that isn't already missing, set to missing based on GQ
  to_filter <- snptable[,paste0(pool,".GQ")] < GQ_threshold
  to_filter[is.na(to_filter)] <- FALSE
  snptable[to_filter, paste0(pool,".GT")] <- "./."
  snptable[to_filter, paste0(pool,c(".DP",".FREQ",".PVAL",".AD",".RD"))] <- NA
  
}

##### Now remove any rows with > 25% missing data... ##### 
cat(paste0(">>> Removing sites with missing data >",max_missing_threshold,"\n"))

missing_count <- rowSums(snptable[,paste0(pool_names,".GT")] == "./.")
missing_to_keep <- missing_count/length(pool_names) < max_missing_threshold
snptable_nomiss <- snptable[missing_to_keep,]

snps_removed <- nrow(snptable) - nrow(snptable_nomiss)
message(paste0(">>> Removed ",snps_removed," of ",nrow(snptable)," SNPs based on missing data of ",max_missing_threshold*100,"%
               "))


##### Now calculate the global allele frequency #####
# Fetch the ploidy information
cat(paste0(">>> Removing sites with minor allele frequency < ",maf_threshold,"\n"))
# if(is.null(ploidy_path)){
#   ploidy_file <- data.frame(pool=pool_names,
#                             ploidy=1)
# } else {

# Take the weighted average across the pools, weighted by ploidy
all_freqs <- snptable_nomiss[,paste0(pool_names,".FREQ")]
all_freqs <- apply(all_freqs,2,gsub,pattern="%",replacement="")
all_freqs <- apply(all_freqs,2,gsub,pattern=",",replacement=".")
class(all_freqs) <- "numeric"

# Calculate the global frequency of the REF allele
global_freqs <- (apply(all_freqs,1,weighted.mean,w=ploidy_file$ploidy,na.rm=T))/100

# Add this column to the snptable
snptable_nomiss$locus <- paste0(snptable_nomiss$CHROM,"-",snptable_nomiss$POS)
global_maf <- global_freqs
global_maf[global_maf > 0.5] <- 1 - global_maf[global_maf > 0.5]
snptable_nomiss$MAF <- global_maf

# Filter for maf
snptable_nomiss_nomaf <- snptable_nomiss[snptable_nomiss$MAF > maf_threshold,]
snps_removed <- nrow(snptable_nomiss) - nrow(snptable_nomiss_nomaf)
message(paste0(">>> Removed ",snps_removed," of ",nrow(snptable_nomiss)," SNPs based on MAF of ",maf_threshold,"
               "))

# Report final snps
message(paste0(">>> Final SNP Count = ",nrow(snptable_nomiss_nomaf)-1,"
               "))

# And write this to an output
output_file <- gsub("unfiltered",paste0("filtered_NoMiss",max_missing_threshold,"_maf",maf_threshold),snptable_path)
cat(paste0(">>> Finished, writing filtered output to ",output_file,"\n"))
write.table(snptable_nomiss_nomaf,
            output_file,
            row.names = F,quote = F,sep = "\t")

