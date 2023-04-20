####################################################################################
# This script takes a VCF input and performs analyses of pop structure: 
# Summarise structure as: PCA + Global FST + Spatial Structure

# Load these
lib <- c("vegan","geosphere","R.utils","tidyverse","doParallel","poolfstat","parallel","data.table","sp","VGAM","argparse","vcfR","gdsfmt","SNPRelate","ggplot2")
dummy <- suppressPackageStartupMessages(lapply(lib,function(x){ 
  if (!require(x,character.only = T)){ 
    install.packages(x,repos='https://utstat.toronto.edu/cran/')
    library(x,character.only = T)
  } else {
    library(x,character.only = T)
  }
}))

# Set a seed
set.seed(1000)

##########################################################################################
# # Set up our arguments to parse
# Read the VCF in from the command line
args <- commandArgs(asValues = T,excludeReserved = T)[-1]
print(args)

vcf_path <- args$vcf
n_cores <- as.integer(args$n_cores)
snp_downsample <- as.integer(args$sub_SNP)
is.pool <- args$pool
output_dir <- args$dataset_dir


if(any(is.null(c(vcf_path,output_dir)))){
  stop("ERROR: No VCF/OUTPUT provided")
}
if(is.null(n_cores)){
  n_cores <- 1
}
if(is.null(snp_downsample)){
  snp_downsample <- 10000
}
if(is.null(is.pool)){
  is.pool <- FALSE
}

# #### For debugging...
# vcf_path <- "data/VCFs/14b_Athaliana_Weigel_wgs_IBE/weigel_Athaliana_IBE_full_concatened_maf05.vcf.gz"
# n_cores <- 16
# is.pool <- FALSE
# output_dir <- "test_new_pipeline_weigel_IBE"
# snp_downsample=10000

# Set up parallel cores
cl <- makeCluster(n_cores, type="FORK")
registerDoParallel(cl)

# Kill cluster
on.exit(stopCluster(cl))

################################################################################################
# Based on the VCF, read in the metadata and subset it if needs be...
if(file.exists(paste0(dirname(vcf_path),"/sampling_data.csv"))){
  metadata <- read.csv(paste0(dirname(vcf_path),"/sampling_data.csv"))
} else {
  stop("Error: Expecting a file with sampling info as 'sampling_data.csv' in the same directory as the VCF, but this does not exist...")
}

# Make if not made for whatever reason...
dir.create("outputs/GEA_res",showWarnings = F)

# Fetch VCF individuals and further subset...
if(is.pool){
  
  # Pull column headers and clean
  vcf_inds <- system(paste0("head -n1 ",vcf_path),intern = T)
  vcf_inds <- unlist(strsplit(vcf_inds,"\t"))
  vcf_inds <- gsub(".GT","",grep(".GT",vcf_inds,value = T))
  
} else {
  vcf_inds <- system(paste0("bin/bcftools query -l ",vcf_path),intern = T)
}

# Do we need to clean?
if(any(!(vcf_inds %in% metadata$Sample_name))){
  to_remove <- c(vcf_inds,metadata$Sample_name)[!(duplicated(c(vcf_inds,metadata$Sample_name)))]
  message(paste0(">>> Warning: The following individuals/pools are missing from VCF or sampling data: ",paste(to_remove,collapse=",")))
  message(paste0(">>> Warning: No additional SNP filtering will be done to account for missing individuals"))
  
  metadata <- metadata[!(metadata$Sample_name %in% to_remove),]
  vcf_inds <- vcf_inds[!(vcf_inds %in% to_remove)]
  
}

# List final vcf inds
vcf_inds_final <- vcf_inds

# Make directory and set
species_code <- output_dir
workdir <- paste0("outputs/GEA_res/",species_code)
suppressWarnings(dir.create(workdir))

# Set up geo pops
metadata$geo_pops <- paste0(metadata$Lat,"_",metadata$Long)


# Now calculate population structure --------------------------------------
# IF the data is pooled, check that we don't have duplicated geopops and if we do, bump them...
if(is.pool & any(duplicated(metadata$geo_pops))){
  message(">>> WARNING: Duplicated geopops detected. Only taking the first population for each geo-coord.
          ")
  to_remove <- duplicated(metadata$geo_pops)
  metadata <- metadata[!(to_remove),]
  
  # List final vcf inds
  vcf_inds_final <- vcf_inds[vcf_inds %in% metadata$Sample_name]
}

# Report
message(paste0(">>> VCF contains: 
               Species = ",species_code,"
               Inds = ",length(vcf_inds_final),"
               Pops = ",length(unique(metadata$geo_pop)),"
               "))


if(!is.pool){
  message(">> VCF is WGS/Capture
          ")
  message(">>> Parsing full VCF into SNPRelate
          ")
  # Make a temp file for snpgds conversion...
  snpgds_tempfile <- tempfile(pattern = "snpgds_tmp", fileext = '.gds')
  on.exit({ unlink(snpgds_tempfile) })
  
  # Convert
  snpgdsVCF2GDS(vcf_path,snpgds_tempfile,method = "biallelic.only")
  
  # Read in
  gds_input <- snpgdsOpen(snpgds_tempfile)
  
  # Linkage prune using our current individuals
  message(">>> Pruning for linkage of r2 0.2 with SNPRelate with 50kb windows
          ")
  snpset <- snpgdsLDpruning(gds_input, 
                            sample.id = vcf_inds_final, 
                            ld.threshold=0.2,
                            autosome.only = F,
                            num.thread = n_cores,
                            slide.max.bp=50000L)
  
  # Subset for our unlinked SNPs
  if(length(unlist(snpset)) > snp_downsample){
    random_unlinked <- sample(unlist(snpset),snp_downsample)
  } else {
    random_unlinked <- unlist(snpset)
  }
  
  # Perform PCA
  message(">>> Performing PCA and Fst structure assessment
          ")
  pca <- snpgdsPCA(gds_input, sample.id = vcf_inds_final, snp.id=random_unlinked, num.thread=n_cores,autosome.only = F,eigen.cnt = 0)
  
  # Primary eigenvector
  eig1 <- pca$varprop[1]
  
  # 50% variance
  eig_cumsum <- cumsum(pca$varprop)
  eig50 <- length(eig_cumsum[eig_cumsum < 0.5]) + 1
  
  # Now do FST based on population assignments...
  popmap <- metadata[,c("Sample_name","geo_pops")]
  snprelate_inds <- read.gdsn(index.gdsn(gds_input, "sample.id"))
  snprelate_inds <- snprelate_inds[snprelate_inds %in% popmap$Sample_name]
  popmap <- popmap[match(snprelate_inds,popmap$Sample_name),]
  
  # Produce a figure of PCA
  pca_res <- data.frame(ind=pca$sample.id,
                        PC1=pca$eigenvect[,1],
                        PC2=pca$eigenvect[,2])
  pca_res$pop <- popmap$geo_pops
  pca_fig <- ggplot(pca_res,aes(PC1,PC2,colour=pop))+
    geom_point()+
    theme_minimal()+
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=16),
          legend.position = "none")+
    labs(x=paste0("PC1 (",round(pca$varprop[1]*100,2),"%)"),
         y=paste0("PC2 (",round(pca$varprop[2]*100,2),"%)"))+
    ggtitle(gsub("outputs/","",species_code))
  
  # Save to results dir
  pdf(paste0(workdir,"/SNPRelate_pca_fig.pdf"),width=10,height=8)
  print(pca_fig)
  dev.off()
  
  # Caculate W+C 1984 Fst
  species_fst <- snpgdsFst(gds_input,
                           sample.id = popmap$Sample_name,
                           snp.id = random_unlinked, 
                           population=as.factor(popmap$geo_pops),
                           method="W&C84",autosome.only = F)
  
  # Caculate W+H 2002 Fst for beta metrix
  pop_fst <- snpgdsFst(gds_input,
                       sample.id = popmap$Sample_name,
                       snp.id = random_unlinked, 
                       population=as.factor(popmap$geo_pops),
                       method="W&H02",autosome.only = F)
  
  # Remove diagonal
  pop_fst$Beta[diag(pop_fst$Beta)] <- NA
  
  # We also want to calculate the physical distance between populations...
  all_pops <- sort(unique(popmap$geo_pops))
  pop_locations <- data.frame(all_pops=all_pops) %>% tidyr::separate("all_pops",into=c("Lat","Long"),sep="_")
  pop_locations <- as.matrix(pop_locations[,c("Long","Lat")])
  class(pop_locations) <- "numeric"
  pop_dist_mat <- distm(pop_locations, fun = distHaversine)
  colnames(pop_dist_mat) <- rownames(pop_dist_mat) <- all_pops
  
  # Mantel test of Fst vs Distance
  fst_dist_mantel <- mantel(pop_fst$Beta/(1-pop_fst$Beta), pop_dist_mat, method = "kendall", permutations = 9999, na.rm = TRUE)
  
  # Save everything as a list and an RDS
  saveRDS(list(pca=pca,
               fst=species_fst,
               eig1=eig1,
               eig50=eig50,
               fst_dist=fst_dist_mantel,
               fst_mat=pop_fst$Beta,
               pop_dist_mat=pop_dist_mat,
               fst_dist_mantel=fst_dist_mantel),
          paste0(workdir,"/SNPRelate_pca_fst_results.rds"))
  
} else {
  message(">>> VCF is PoolSeq, using AFs
          ")
  
  # Read in
  pool_AFs <- data.frame(fread(vcf_path,select=c("CHROM","POS",paste0(vcf_inds_final,".FREQ")),header=T))
  colnames(pool_AFs) <- c("chr","pos",vcf_inds_final)
  
  # Clean away %
  pool_AFs[,vcf_inds_final] <- apply(pool_AFs[,vcf_inds_final],2,parse_number)
  for(ind in vcf_inds_final){
    pool_AFs[,ind] <- pool_AFs[,ind]/100
  }
  
  # Sample
  pool_AFs_sub <- pool_AFs[sort(sample(1:nrow(pool_AFs),snp_downsample)),]
  
  # Build corr mat and check for no correlations
  LD_cor <- cor(t(pool_AFs_sub[,vcf_inds_final]))^2
  diag(LD_cor) <- NA
  
  # Record the proportion of remaining SNPs over LD R^2 of 0.2...
  ld_vector <- na.omit(as.numeric(LD_cor))
  linkage_prop <- length(ld_vector[ld_vector > 0.2])/length(ld_vector)
  
  # Report
  message(paste0(">>> After LD filtering, poolseq VCF has ",round(linkage_prop*100,3),"% of SNPs in LD (r2) > 0.2
                 "))
  
  # Now perform PCA based on frequencies for subset
  # Replace any missing values
  for(i in 1:nrow(pool_AFs_sub)){
    if(any(is.na(pool_AFs_sub[i,vcf_inds_final]))){
      
      to_change <- is.na(pool_AFs_sub[i,vcf_inds_final])
      to_change <- colnames(to_change)[to_change]
      pool_AFs_sub[i,to_change] <- rowMeans(pool_AFs_sub[i,vcf_inds_final],na.rm = T)
    }
  }
  pool_pca <- prcomp(t(pool_AFs_sub[,vcf_inds_final]),center = T,scale. = F)
  
  # Collect summaries of PCA
  # Primary eigenvector
  eigenvals <- pool_pca$sdev^2/sum(pool_pca$sdev^2)
  eig1 <- eigenvals[1]
  
  # 50% variance
  eig_cumsum <- cumsum(eigenvals)
  eig50 <- length(eig_cumsum[eig_cumsum < 0.5]) + 1
  
  # Calculate global Fst from a subsetted VCF
  SNPs_to_keep <- pool_AFs_sub[,1:2]
  write.table(SNPs_to_keep,paste0(workdir,"/subset_SNPs_LD.txt"),row.names = F,col.names = F,sep = "\t",quote = F)
  
  # Make a baypass counts file
  # REF DEPTH
  pool_RD <- data.frame(fread(vcf_path,select=c("CHROM","POS",paste0(vcf_inds_final,".RD")),header=T))
  pool_RD <- pool_RD[paste0(pool_RD$CHROM,"_",pool_RD$POS) %in% paste0(SNPs_to_keep[,1],"_",SNPs_to_keep[,2]),]
  
  # ALT DEPTH
  pool_AD <- data.frame(fread(vcf_path,select=c("CHROM","POS",paste0(vcf_inds_final,".AD")),header=T))
  pool_AD <- pool_AD[paste0(pool_AD$CHROM,"_",pool_AD$POS) %in% paste0(SNPs_to_keep[,1],"_",SNPs_to_keep[,2]),]
  
  # Make a baypass format allele count file
  baypass_out <- matrix(ncol=2*length(vcf_inds_final),nrow=nrow(pool_AD))
  oddN <- seq(1,ncol(baypass_out),2)
  evenN <- seq(2,ncol(baypass_out),2)
  for(i in 1:length(oddN)){
    baypass_out[,oddN[i]] <- pool_RD[,paste0(vcf_inds_final[i],".RD")]
    baypass_out[,evenN[i]] <- pool_AD[,paste0(vcf_inds_final[i],".AD")]
  }
  
  # Write this out to the workdir with missing data coded according to manual...
  baypass_out[is.na(baypass_out)] <- 0
  write.table(baypass_out,
              paste0(workdir,"/snp_downsample_baypass_format.txt"),
              sep="\t",col.names = F,row.names = F,quote = F)
  
  # Also make a haploid size file from meta
  all_haploids <- metadata[metadata$Sample_name %in% vcf_inds_final,c("Sample_name","Pool_Size")]
  all_haploids <- all_haploids[match(all_haploids$Sample_name,vcf_inds_final),]
  haploid_size_mat <- matrix(ncol=length(vcf_inds_final),nrow=1)
  haploid_size_mat[1,] <- all_haploids$Pool_Size
  
  # If poolsize is empty, assume pools of 20 individuals...
  if(all(is.na(haploid_size_mat[1,]))){
    haploid_size_mat[1,] <- 20
  }
  
  write.table(haploid_size_mat,
              paste0(workdir,"/haploid_size_baypass_format.txt"),
              sep="\t",col.names = F,row.names = F,quote = F)
  
  # Now read back in poolformat
  poolvcf <- genobaypass2pooldata(
    genobaypass.file = paste0(workdir,"/snp_downsample_baypass_format.txt"),
    poolsize.file = paste0(workdir,"/haploid_size_baypass_format.txt"),
    snp.pos = SNPs_to_keep,
    poolnames = vcf_inds_final,
    min.cov.per.pool = -1,
    max.cov.per.pool = 1e+06,
    verbose = TRUE
  )
  
  # Get global Fst
  species_fst <- computeFST(poolvcf)
  
  # Report
  message(">>> Pool Global Fst of ",round(species_fst$FST,3),"
          ")
  
  # Get pairwise fst
  pop_fst <- compute.pairwiseFST(poolvcf,verbose=FALSE)
  
  # Fix row names...
  for(i in 1:ncol(pop_fst@PairwiseFSTmatrix)){
    rownames(pop_fst@PairwiseFSTmatrix)[i] <- colnames(pop_fst@PairwiseFSTmatrix)[i] <- paste(metadata[metadata$Sample_name==colnames(pop_fst@PairwiseFSTmatrix)[i],c("Lat","Long")],collapse = "_")
  }
  
  
  # We also want to calculate the physical distance between populations...
  all_pops <- sort(paste0(metadata$Lat,"_",metadata$Long))
  pop_locations <- data.frame(all_pops) %>% separate("all_pops",into=c("Lat","Long"),sep="_")
  pop_locations <- as.matrix(pop_locations[,c("Long","Lat")])
  class(pop_locations) <- "numeric"
  pop_dist_mat <- distm(pop_locations, fun = distHaversine)
  colnames(pop_dist_mat) <- rownames(pop_dist_mat) <- all_pops
  
  # Fix order...
  pop_fst@PairwiseFSTmatrix <- pop_fst@PairwiseFSTmatrix[colnames(pop_dist_mat),colnames(pop_dist_mat)]
  
  # Mantel test of Fst vs Distance
  fst_dist_mantel <- mantel(pop_fst@PairwiseFSTmatrix/(1-pop_fst@PairwiseFSTmatrix), pop_dist_mat, method = "kendall", permutations = 9999, na.rm = TRUE)
  
  # Save the results
  saveRDS(list(pca=pool_pca,
               fst=species_fst,
               eig1=eig1,
               eig50=eig50,
               linkage_prop=linkage_prop,
               fst_dist=fst_dist_mantel,
               fst_mat=pop_fst@PairwiseFSTmatrix,
               pop_dist_mat=pop_dist_mat,
               fst_dist_mantel=fst_dist_mantel),
          paste0(workdir,"/SNPRelate_pca_fst_results.rds"))
  
  # Tidy
  system(paste0("rm -f ",workdir,"/subsampled_poolvcf.recode.vcf"))
}

# Send completion
message(">>> All Finished!
        ")




