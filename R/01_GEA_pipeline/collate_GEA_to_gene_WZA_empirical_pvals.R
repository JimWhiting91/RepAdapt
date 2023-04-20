####################################################################################
# This script takes our GEA results and condenses them down into per gene WZA and Top-Candidate (TC) scores
lib <- c("pbmcapply","doParallel","R.utils","argparse","parallel","data.table","sp","VGAM","tidyr","dplyr")
dummy <- suppressPackageStartupMessages(lapply(lib,function(x){ 
  if (!require(x,character.only = T)){ 
    install.packages(x,repos='https://cloud.r-project.org/')
  } else {
    library(x,character.only = T)
  }
}))

# Read the VCF in from the command line
args <- commandArgs(asValues = T,excludeReserved = T)[-1]
print(args)

vcf_path <- args$vcf
gff_path <- args$gff
n_cores <- as.integer(args$n_cores)
TC_threshold <- as.numeric(args$tc_threshold)
output_dir <- args$dataset_dir
snp_per_gene <- as.numeric(args$snp_per_gene)
gene_flank_size <- as.integer(args$gene_flank_size)

if(any(is.null(c(vcf_path,gff_path,output_dir)))){
  stop("ERROR: No VCF/GFF/OUTPUT provided")
}
if(is.null(n_cores)){
  n_cores <- 1
}
if(is.null(TC_threshold)){
  TC_threshold <- 0.99
}
if(is.null(snp_per_gene)){
  snp_per_gene <- 20
}
if(is.null(gene_flank_size)){
  gene_flank_size <- 500
}

# ################################################
# # For setting up purposes
# vcf_path <- "data/VCFs/22_Pabies_Milesi/milesi_Pabies_full_concatened_maf05_unstitched.vcf.gz"
# gff_path <- "data/reference_genomes/P.abies_reference/Pabies01-gene_CHR_MATCHED.gff3"
# TC_threshold <- 0.99
# n_cores <- 8
# output_dir <- "230301_Picea_abies_Milesi_Capture"
# gene_flank_size <- 500
# snp_per_gene = 0.75

# Set up workdir...
workdir <- paste0("outputs/GEA_res/",output_dir)
suppressWarnings(dir.create("outputs/GEA_res"))
suppressWarnings(dir.create(workdir))

climate_vars = gsub("_GEA.rds","",list.files(workdir,pattern = "_GEA.rds"))

################################################
##### PREPARE GENE SNP OVERLAPS #####
message(">>> Building keys for SNPs > genes
        ")

# Read in the start-end positions based on the file linking the GFF and the OF from the gff directory...
gff <- read.table(paste0(dirname(gff_path),"/",list.files(dirname(gff_path),pattern = "_proteome_to_OF_id_map.txt")),header=T,sep="\t")
colnames(gff) <- c("gene_ID","OF_ID","seqname","start","end","gea_gene")
gff$gea_gene = paste0(gff$seqname,":",gff$start,"-",gff$end)
# Set up bed file for all genes
gene_dd <- data.table(chr=gff$seqname,
                      start=gff$start,
                      end=gff$end,
                      OF_ID = gff$OF_ID)

# Add a 500bp flank to all genes...
gene_dd$start <- ifelse(gene_dd$start - gene_flank_size < 0,0,gene_dd$start - gene_flank_size)
gene_dd$end <- gene_dd$end + gene_flank_size


# Set up keyed index for overlap
gene_bed <- gene_dd[,.(chr,start,end,OF_ID)]
setkey(gene_bed, chr, start, end)

# Set up SNP boundaries
GEA_res <- list.files(workdir,pattern="GEA.rds")
if(grepl(".vcf",vcf_path)){
  GEA_SNPs <- system(paste0("bcftools query -f '%CHROM\t%POS\n' ",vcf_path),intern = T)
} else {
  GEA_SNPs <- system(paste0("tail -n+2 ",vcf_path," | cut -f1,2"),intern = T)
}

# Fetch chromosomes and SNPs
snp_dd <- data.frame(chr=sapply(strsplit(GEA_SNPs,"\t"),'[[',1),
                     start=as.integer(sapply(strsplit(GEA_SNPs,"\t"),'[[',2)))

# Set these up as genomic regions
snp_dd$end <- snp_dd$start+1
snp_bed <- data.table(snp_dd)
setkey(snp_bed, chr, start, end)

# Find overlaps
snp_gene_overlap <- foverlaps(snp_bed, gene_bed, minoverlap=1L,nomatch=NA)

# Use overlaps to group genes
snp_gene_overlap$snp_id <- paste0(snp_gene_overlap$chr,":",snp_gene_overlap$i.start)
snp_gene_overlap$gene_id <- paste0(snp_gene_overlap$chr,":",snp_gene_overlap$start,"-",snp_gene_overlap$end)
snp_gene_overlap[is.na(snp_gene_overlap$start),gene_id := NA]

##########################################################################################
# Calculate pbar-qbar -----------------------------------------------------
# Read in all the frequency files to fetch pbar-qbar
# Frequency files are calculated prior to GEA in perform_GEA.R
message(">>> Calculating pbar-qbar
        ")
af_res <- readRDS(paste0(workdir,"/pop_allele_frqs.rds"))

# Build up p_mat and q_mat
p_mat <- matrix(nrow=nrow(af_res[[1]]),ncol=length(af_res))
q_mat <- p_mat
for(i in 1:ncol(p_mat)){
  p_mat[,i] <- af_res[[i]]$ref_freq
  q_mat[,i] <- af_res[[i]]$alt_freq
}

# And calculate pbar-qbar
pbar_qbar <- data.table(snp_id=paste0(af_res[[1]]$chr,":",af_res[[1]]$bp),
                        pbar_qbar = rowMeans(p_mat,na.rm = T) * rowMeans(q_mat,na.rm = T))

################################################
##### ESTIMATE WZA FOR EACH GENE #####
# Announce
message(">>> Calculating WZA with empirical Pvals for all climate variables
        ")
for(res in GEA_res){

  #########################################################################################################################
  #### WZA approach #####
  
  if(!file.exists(paste0(workdir,"/",gsub("_GEA.rds","_WZA_TC_allgenes.rds",res)))){
    
    # Fetch res
    GEA <- data.table(readRDS(paste0(workdir,"/",res)))
    
    # Remove the non-snp genes...
    snps_to_keep = snp_gene_overlap[!is.na(gene_id),snp_id]
    
    # Add pbar-qbar 
    GEA_snps <- merge(GEA[snp_id %in% snps_to_keep,],
                 pbar_qbar[snp_id %in% snps_to_keep,],by="snp_id")
    
    # Merge GEA and snp_gene_overlap
    GEA_merge <- merge(snp_gene_overlap[snp_id %in% snps_to_keep,],
                       GEA_snps,by="snp_id")
    
    # Filter away SNPs with NA pvals, or pvals=1
    GEA_merge <- GEA_merge[!(is.na(GEA_merge$p)),]
    
    # Convert to empirical Pvals
    # Adding the +1 here just to catch any dodgy cases where there is 1/1
    GEA_merge$empirical_p <- rank(GEA_merge$p)/(nrow(GEA_merge)+1)
    
    # Convert one-sided p-values to Z-scores
    GEA_merge$z <- qnorm(GEA_merge$empirical_p , lower.tail = F)
    
    ## Calculate the numerator of the Weighted-Z score
    weiZ_num <- tapply(GEA_merge$pbar_qbar * GEA_merge$z, GEA_merge$gene_id, sum )
    
    ## Calculate the denominator of the Weighted-Z score
    weiZ_den <- sqrt(tapply( GEA_merge$pbar_qbar^2, GEA_merge$gene_id, sum ))
    
    ## Bring data together into a new DF
    Z_df <- data.table(gene_id = names(weiZ_num), weiZ = weiZ_num/weiZ_den)
    
    #########################################################################################################################
    #### Repeat WZA with downsampling #####
    
    # Re-estimate snp_per_gene if we are using a proportion
    if(snp_per_gene < 1){
      snp_per_gene <- quantile(table(GEA_merge$gene_id),probs=snp_per_gene)
      
      # But make sure that we cap this at a threshold of 10
      if(snp_per_gene < 10){
        snp_per_gene <- 10
      }
    }
    
    # For each gene, add a random 1:nrow and then only retain the first N
    message(paste0(">>> Starting WZA downsample for ",GEA$climate_var[1]," for ",workdir,"
                 "))
    gene_counts = table(GEA_merge$gene_id)
    genes_to_downsample = names(gene_counts)[gene_counts >= snp_per_gene]
    GEA_merge_downsample = GEA_merge[gene_id %in% genes_to_downsample,]
    wza_downsample <- rbindlist(lapply(1:100,function(iter){
      set.seed(iter)
      
      # Randomly choose up to snp_per_gene SNPs per gene
      GEA_merge_sample <- GEA_merge_downsample[sample(1:nrow(GEA_merge_downsample)),]
      GEA_merge_sample <- GEA_merge_sample[,.(snp_id,
                                              pbar_qbar,
                                              z,
                                              snp_i=1:nrow(.SD)),by=gene_id]
      GEA_merge_sample <- GEA_merge_sample[GEA_merge_sample$snp_i <= snp_per_gene,]
      
      # Calculate WZA as before...
      ## Calculate the numerator of the Weighted-Z score
      weiZ_num_sample <- tapply(GEA_merge_sample$pbar_qbar * GEA_merge_sample$z, GEA_merge_sample$gene_id, sum )
      
      ## Calculate the denominator of the Weighted-Z score
      weiZ_den_sample <- sqrt(tapply( GEA_merge_sample$pbar_qbar^2, GEA_merge_sample$gene_id, sum ))
      
      ## Bring data together into a new DF
      Z_df_tmp <- data.table(gene_id = names(weiZ_num_sample),
                         weiZ = weiZ_num_sample/weiZ_den_sample,
                         iteration=iter)
      
    }))
    
    wza_downsample_final <- wza_downsample[,.(weiZ_downsample=mean(.SD$weiZ)),by=gene_id]
    
    # Merge this with the full WZA...
    Z_df_merge <- merge(Z_df,wza_downsample_final,by="gene_id",all.x = T)
    Z_df_merge[is.na(weiZ_downsample),weiZ_downsample := weiZ]
    
    #########################################################################################################################
    #### TC approach #####
    
    # Announce
    message(paste0(">>> Calculating TC for ",GEA$climate_var[1]," for ",workdir,"
                 "))
    
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
    
    # Finally merge with weiZ, format, and save
    final_res_merge <- merge(Z_df_merge,TC_merge_final,by="gene_id")

    # Also fetch the mean pbar_qbar
    mean_pbar_qbar <- suppressMessages(
      data.frame(GEA_merge %>% group_by(gene_id) %>% summarise(mean_pbar_qbar=mean(pbar_qbar)))
    )
    final_res_merge <- merge(final_res_merge,mean_pbar_qbar)
    
    # Make sure we know how many SNPs were used for downsampling
    final_res_merge$downsampled_snp_count <- final_res_merge$snp_count
    final_res_merge$downsampled_snp_count[as.integer(final_res_merge$downsampled_snp_count) > as.integer(snp_per_gene)] <- snp_per_gene
    
    # And merge with original OF IDs
    final_res_merge_OF = merge(final_res_merge,unique(GEA_merge[,.(gene_id,OF_ID)]),by = "gene_id")
    
    # Save the results
    saveRDS(final_res_merge_OF,
            paste0(workdir,"/",unique(final_res_merge$climate_var),"_WZA_TC_allgenes.rds"))
  }
}

# Send completion
message(">>> All Finished!
        ")




