####################################################################################
# This script takes a VCF input and performs GEA analysis over climate variables
lib <- c("parallel","data.table","sp","VGAM","tidyr","dplyr")
sapply(lib,library,character.only=T)

# Read the VCF in from the command line
args <- commandArgs(TRUE)
vcf_path <- as.character(args[1])
metadata_path <- as.character(args[2])
n_cores <- as.integer(args[3])
gff_path <- as.character(args[4])
TC_threshold <- as.numeric(args[5])

# # For setting up purposes
# vcf_path <- "/lu213/james.whiting/RepAdapt/data/VCFs/06_Ahalleri_Kubota_wgs/Ahalg_full_concatened.vcf.gz"
# results_dir <- "outputs/Arabidopsis_halleri_Kubota_Individual/"
# gff_path <- "data/reference_genomes/A.halleri_gemmifera_reference/Arabidopsis_halleri.Ahal2.2.45.gff3"
# metadata_path <- "metadata/sample_species_vcf_author_map_v2_210507.csv"
# TC_threshold <- 0.99
# 
# n_cores <- 19

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
snp_gene_overlap <- data.frame(foverlaps(snp_bed, gene_bed, minoverlap=1L,nomatch=NA))

# Use overlaps to group genes
snp_gene_overlap$snp_id <- paste0(snp_gene_overlap$chr,":",snp_gene_overlap$i.start)
snp_gene_overlap$gene_id <- paste0(snp_gene_overlap$chr,":",snp_gene_overlap$start,"-",snp_gene_overlap$end)
snp_gene_overlap[is.na(snp_gene_overlap$start),"gene_id"] <- NA

# Remove away empty SNP rows 
# gene_snp_overlap <- na.omit(gene_snp_overlap)

# # And add these to original beds
# gene_dd$gene_id <- paste0(gene_dd$chr,":",gene_dd$start,"-",gene_dd$end)
# snp_dd$snp_id <- paste0(snp_dd$chr,":",snp_dd$start)
# 
# # Make a list of SNP IDs for relevant genes...
# gene_snp_list <- pbmcapply::pbmclapply(unique(gene_snp_overlap$gene_id), function(x) return(gene_snp_overlap[gene_snp_overlap$gene_id == x,"snp_id"]),mc.cores=n_cores)

################################################
##### LOOP OVER DATASETS WITHIN THE VCF #####

# Fetch info on which species we are looking at and set working directory
metadata <- read.csv(metadata_path)
sub_meta <- metadata[metadata$VCF==vcf_path,]

# Make a working dir for each species that shares the VCF
sub_meta$Taxon_name <- gsub(" ","_",sub_meta$Taxon_name)
sub_meta$species_codes <- paste0("outputs/",sub_meta$Taxon_name,"_",sub_meta$Author,"_",sub_meta$Data_type)

# Fetch VCF individuals and further subset...
vcf_inds <- system(paste0("bin/bcftools query -l ",vcf_path),intern = T)
sub_meta <- sub_meta[sub_meta$Sample_name %in% vcf_inds,]

# Remove species that appear <10 times...
species_counts <- table(sub_meta$species_codes)

# Report on this
species_codes <- names(species_counts[species_counts >= 10])

# Loop over species code and summarise GEA
for(results_dir in species_codes){
  
  ################################################
  ##### PREPARE ALLELE FREQUENCY INPUT #####
  
  # # Based on the VCF, read in the metadata and subset it
  # metadata <- read.csv(metadata_path)
  # sub_meta <- metadata[metadata$VCF==vcf_path,]
  # 
  # # Get the individuals we're getting frequencies for
  # inds <- sub_meta[,"Sample_name"]
  # inds_to_keep <- paste0("--indv ",inds,collapse = " ")
  # 
  # # Calculate global allele frequencies with vcftools only keeping individuals with metadata
  # system(paste0("vcftools --gzvcf ",vcf_path," --freq ",inds_to_keep," --out ",results_dir,"/whole_data_AFs"),wait = T)
  
  # Read in all the frequency files to fetch pbar-qbar
  freq_files <- list.files(results_dir,pattern=".frq")
  
  # Fetch all the minor and major allele frequencies
  af_res <- mclapply(freq_files,function(freq){
    
    # Fetch the allele frequency results
    freq_res <- data.frame(fread(paste0(results_dir,"/",freq),fill=T))
    colnames(freq_res) <- c("chr","bp","allele_N","allele_count","ref_freq","alt_freq")
    
    # Format for analysis
    freq_dd <- data.frame(chr=freq_res$chr,
                          bp=freq_res$bp,
                          ref_freq=as.numeric(sapply(strsplit(freq_res$ref_freq,":"),'[[',2)),
                          alt_freq=as.numeric(sapply(strsplit(freq_res$alt_freq,":"),'[[',2)))
    
    # And make a vector of minor allele frequencies
    # freq_dd$minor_af <- rje::rowMins(freq_dd[,c("ref_freq","alt_freq")])
    return(freq_dd)
  },mc.cores=n_cores)
  
  # Build up p_mat and q_mat
  p_mat <- matrix(nrow=nrow(af_res[[1]]),ncol=length(af_res))
  q_mat <- p_mat
  for(i in 1:ncol(p_mat)){
    p_mat[,i] <- af_res[[i]]$ref_freq
    q_mat[,i] <- af_res[[i]]$alt_freq
  }
  
  # And calculate pbar-qbar
  pbar_qbar <- data.frame(snp_id=paste0(af_res[[1]]$chr,":",af_res[[1]]$bp),
                          pbar_qbar = rowMeans(p_mat,na.rm = T) * rowMeans(q_mat,na.rm = T))
  
  # # Remove species identifier if it is present
  # pbar_qbar$snp_id <- sapply(strsplit(pbar_qbar$snp_id,"_"),'[[',2)
  
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
  
  # Run for all GEA outputs
  GEA_res <- list.files(results_dir,pattern="GEA.tsv")
  
  # Apply over all GEA
  mclapply(GEA_res,function(res){
    
    #########################################################################################################################
    #### WZA approach #####
    
    # Fetch res
    GEA <- data.frame(fread(paste0(results_dir,"/",res)))
    colnames(GEA)[colnames(GEA) == "snp"] <- "snp_id"
    
    # Announce
    message(paste0("Calculating WZA for ",GEA$climate_var[1]," for ",results_dir))
    
    # Clean SNP col
    GEA$snp_id <- sapply(strsplit(GEA$snp_id,"_"),'[[',2)
    
    # Add pbar-qbar 
    GEA$pbar_qbar <- pbar_qbar$pbar_qbar
    
    # Merge GEA and snp_gene_overlap
    GEA_merge <- merge(snp_gene_overlap,GEA,by="snp_id")
    
    # Return only SNPs located in genes
    GEA_merge <- GEA_merge[!(is.na(GEA_merge$gene_id)),]
    
    # Filter away SNPs with NA pvals, or pvals=1
    GEA_merge <- GEA_merge[!(is.na(GEA_merge$p)),]
    GEA_merge[GEA_merge$p == 1, "p"] <- ifelse(max(GEA_merge[GEA_merge$p != 1,"p"]) < 0.99 ,0.99,max(GEA_merge[GEA_merge$p != 1,"p"]))
    
    # Convert one-sided p-values to Z-scores
    GEA_merge$z <- qnorm(GEA_merge$p , lower.tail = F)
    
    ## Calculate the numerator of the Weighted-Z score
    weiZ_num <- tapply(GEA_merge$pbar_qbar * GEA_merge$z, GEA_merge$gene_id, sum )
    
    ## Calculate the denominator of the Weighted-Z score
    weiZ_den <- sqrt(tapply( GEA_merge$pbar_qbar^2, GEA_merge$gene_id, sum ))
    
    ## Bring data together into a new DF
    Z_df <- data.frame(gene_id = names(weiZ_num), weiZ = weiZ_num/weiZ_den)
    
    ## Label the DF with climate variable
    Z_df$climate_var <- unique(GEA_merge$climate_var)
    
    # Save the results
    write.table(Z_df,
                paste0(results_dir,"/",unique(GEA_merge$climate_var),"_WZA_pergene.tsv"),
                row.names = F,quote = F,sep = "\t")
    
    #########################################################################################################################
    #### TC approach #####
    
    # Announce
    message(paste0("Calculating TC for ",GEA$climate_var[1]," for ",results_dir))
    
    # Now we calculate Top Candidate outliers
    # For all genes, we want the number of outliers above the 
    TC_genes <- data.frame(table(GEA_merge$gene_id))
    colnames(TC_genes)[1] <- "gene_id"
    
    # Count outlier
    outlier_counts <- data.frame(table(GEA_merge[abs(GEA_merge$tau.corr) >= quantile(abs(GEA_merge$tau.corr),TC_threshold),"gene_id"]))
    colnames(outlier_counts)[1] <- "gene_id"
    
    # Merge
    TC_merge <- merge(TC_genes,outlier_counts,by="gene_id")
    colnames(TC_merge)[2:3] <- c("snp_count","outlier_count")
    colnames(TC_genes)[2] <- c("snp_count")
    TC_genes$outlier_count <- 0
    
    # Final merge
    TC_merge <- rbind(TC_merge,TC_genes[!(TC_genes$gene_id %in% TC_merge$gene_id),])
    
    # Match up the order of both outputs
    TC_merge <- TC_merge[match(Z_df$gene_id, TC_merge$gene_id),]
    
    # Calculate bionmial expectation expectation
    p <- sum(TC_merge$outlier_count)/sum(TC_merge$snp_count)
    
    # Calculate the expectation
    TC_merge$snp_exp <- qbinom(0.999,TC_merge$snp_count,p)
    
    # Derive outliers
    TC_merge$TC_score <- TC_merge$outlier_count - TC_merge$snp_exp
    
    ## Label the DF with climate variable
    TC_merge$climate_var <- unique(GEA_merge$climate_var)
    
    # Get the mean correlation per gene
    pergene_corr <- suppressMessages(
      data.frame(GEA_merge %>% group_by(gene_id) %>% summarise(mean_corr=mean(abs(tau.corr))))
    )
    TC_merge_final <- merge(TC_merge,pergene_corr,by="gene_id")
    
    # Save the results
    write.table(TC_merge_final,
                paste0(results_dir,"/",unique(TC_merge$climate_var),"_TC_pergene.tsv"),
                row.names = F,quote = F,sep = "\t")
    
    #########################################################################################################################
    ##### Summary of genes ######
    # Here we summarise some info on the genes generally
    # Only do it once
    if(res == GEA_res[1]){
      
      # Get gene count
      gene_summary <- data.frame(gene_id=TC_merge$gene_id,
                                 snp_count=TC_merge$snp_count)
      
      # Get average pbar_qbar for all genes...
      mean_pbar_qbar <- suppressMessages(
        data.frame(GEA_merge %>% group_by(gene_id) %>% summarise(mean_pbar_qbar=mean(pbar_qbar)))
      )
      
      # Merge
      gene_summary <- merge(gene_summary,mean_pbar_qbar,by="gene_id")
      
      # Save this as well
      write.table(gene_summary,
                  paste0(results_dir,"/GEA_snpstats_pergene.tsv"),
                  row.names = F,quote = F,sep = "\t")
      
    }
    
  },mc.cores=n_cores)
  
  # End species code loop
}



