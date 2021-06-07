####################################################################################
# This script takes coadaptree snptables and converts back to VCF format
lib <- c("parallel","data.table","tidyr","dplyr","vcfR","pbmcapply")
dummy <- suppressPackageStartupMessages(lapply(lib,function(x){ 
  if (!require(x,character.only = T)){ 
    install.packages(x,repos='https://cloud.r-project.org/')
  } else {
    library(x,character.only = T)
  }
}))

# # Read the VCF in from the command line
# args <- commandArgs(TRUE)
# snptable <- as.character(args[1])
# file_rows <- as.integer(args[2])
# n_cores <- as.integer(args[3])
# output <- as.character(args[4])

# For setup
snptable <- "~/pooja_orthofinder_info/coadaptree/douglas_fir_coastal/DF_pooled-varscan_all_bedfiles_SNP_FDC.txt"
n_cores <- 16
file_rows=9153701
output <- "douglas_fir_test"

############################################################
snptable_chunk2VCF <- function(snptable_chunk){
  
  # Make fix
  fix <- matrix(nrow=nrow(snptable_chunk),ncol=8)
  fix[,1] <- snptable_chunk[,1]
  fix[,2] <- snptable_chunk[,2]
  fix[,3] <- paste0(fix[,1],"-",fix[,2])
  fix[,4] <- snptable_chunk[,3]
  fix[,5] <- snptable_chunk[,4]
  fix[,6] <- snptable_chunk[,6]
  fix[,7] <- snptable_chunk[,8]
  fix[,8] <- snptable_chunk[,7]
  colnames(fix) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO")
  
  # Get inds
  inds <- grep("GT",colnames(snptable_chunk),value=T)
  inds <- gsub(".GT","",inds)
  
  # Build GT
  gt_mat <- matrix(nrow=nrow(snptable_chunk),ncol=length(inds)+1)
  colnames(gt_mat) <- c("FORMAT",inds)
  colnames(gt_mat) <- c("FORMAT",inds)
  gt_mat[,1] <- paste(snp_details,collapse = ":")
  
  # Fill the mat
  for(ind in inds){
    gt_mat[,ind] <- apply(snptable_chunk[,paste0(ind,".",snp_details)],1,paste,collapse=":")
    gt_mat[,ind] <- gsub(" ","",gt_mat[,ind])
    gt_mat[,ind] <- gsub("e","E",gt_mat[,ind])
  }
  
  return(list(fix=fix,gt_mat=gt_mat))
  
}

# Read in the snptable and handle in chunks
chunksize=10000
file_rows=file_rows-1

# First just capture the header
file_header <- colnames(fread(snptable,nrows=0))

# Now produce some ugly metadata
snp_details <- c("GT","GQ","SDP","DP","RD","AD","FREQ","PVAL")
vcf_meta <- c("##fileformat=VCFv4.1",
              "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",                                                             
              "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">",                                                    
              "##FORMAT=<ID=SDP,Number=1,Type=Integer,Description=\"Raw Read Depth as reported by SAMtools\">",                             
              "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Quality Read Depth of bases with Phred score >= 20\">",                  
              "##FORMAT=<ID=RD,Number=1,Type=Integer,Description=\"Depth of reference-supporting bases (reads1)\">",                        
              "##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Depth of variant-supporting bases (reads2)\">",                          
              "##FORMAT=<ID=FREQ,Number=1,Type=String,Description=\"Variant allele frequency\">",                                           
              "##FORMAT=<ID=PVAL,Number=1,Type=String,Description=\"P-value from Fisher's Exact Test\">"
)

# Now read in all chunks and save them to VCF
dir.create("outputs/snptable_convert")
chunk_seq <- seq(0,file_rows,chunksize)

# Fetch the default data
data(vcfR_example)
vcf_template <- vcf
vcf_template@meta <- vcf_meta

# Loop in parallel
pbmclapply(1:length(chunk_seq),function(x){
  
  # Read in chunk
  snptable_chunk <- data.frame(fread(snptable,nrows = chunksize-1,skip = chunk_seq[x]))
  colnames(snptable_chunk) <- file_header
  
  # Transform
  new_vcf <- snptable_chunk2VCF(snptable_chunk)
  
  # Make new VCF
  out_vcf <- vcf_template
  out_vcf@fix <- new_vcf$fix
  out_vcf@gt <- new_vcf$gt_mat
  
  # Write to the output dir
  write.vcf(out_vcf,paste0("outputs/snptable_convert/",output,"_","chunk",x,".vcf.gz"))
  
  # Index it
  system(paste0("tabix -f -p vcf outputs/snptable_convert/",output,"_","chunk",x,".vcf.gz"))
},mc.cores=n_cores)

# Concatenate them all together
chunk_vcf_paths <- list.files("outputs/snptable_convert/",pattern = output)
chunk_vcf_paths <- sapply(1:length(chunk_vcf_paths),function(x){
  paste0("outputs/snptable_convert/",output,"_","chunk",x,".vcf.gz")
})

# # Write them to a file
# write.table(data.frame(chunk_vcf_paths))

chunk_vcf_merge <- paste(chunk_vcf_paths,collapse = " ")
system(paste0("bcftools concat ",chunk_vcf_merge," -a -Oz -o outputs/snptable_convert/",output,"_merged_snptables.vcf.gz"))
