####################################################################################
# This script takes a VCF input and performs GEA analysis over climate variables
# Load these
lib <- c("parallel","data.table","sp","VGAM","argparse","vcfR","gdsfmt","SNPRelate","ggplot2")
dummy <- suppressPackageStartupMessages(lapply(lib,function(x){ 
  if (!require(x,character.only = T)){ 
    install.packages(x,repos='https://utstat.toronto.edu/cran/')
  } else {
    library(x,character.only = T)
  }
}))

##########################################################################################
# # Set up our arguments to parse
# parser <- ArgumentParser()
# parser$add_argument("-v","--input_vcf", type="character",
#                     help = "Path to gzipped VCF")
# parser$add_argument("-m","--metadata", type="character", 
#                     help = "Path for output rds file")
# parser$add_argument("-nc", "--n_cores", type="integer", default=1, 
#                     help="Number of CPU to use [default %(default)s]",
#                     metavar="number")
# parser$add_argument("-p", "--pool", action="store_false", default=FALSE,
#                     help="Is the VCF based on poolseq? [default %(default)]")
# 
# args <- parser$parse_args()
# ###########################################################################################
# # Get args from command
# vcf_path <- args$input_vcf
# metadata_path <- args$metadata
# n_cores <- args$n_cores

# Read the VCF in from the command line
args <- commandArgs(TRUE)
vcf_path <- as.character(args[1])
metadata_path <- as.character(args[2])
n_cores <- as.integer(args[3])
snp_downsample <- as.integer(args[4])

################################################################################################

# # For setting up purposes
# vcf_path <- "/lu213/james.whiting/RepAdapt/data/VCFs/01_Alyrata_Willi_wgs/Alyrl_full_concatened.vcf.gz"
# metadata_path <- "metadata/sample_species_vcf_author_map_v2_210519.csv"
# n_cores <- 16
# snp_downsample <- 10000

################################################################################################
# Function library

################################################################################################

################################################################################################
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

# Report on this
message(paste0("VCF contains the following species:",paste(species_codes,collapse = ",")))

message("Parsing full VCF into SNPRelate")
# Make a temp file for snpgds conversion...
snpgds_tempfile <- tempfile(pattern = "snpgds_tmp", fileext = '.gds')
on.exit({ unlink(snpgds_tempfile) })

# Convert
snpgdsVCF2GDS(vcf_path,snpgds_tempfile,method = "biallelic.only")

# Read in
gds_input <- snpgdsOpen(snpgds_tempfile)

for(species_code in species_codes){
  #species_code <- species_codes[1]
  
  # Make directory and set
  suppressWarnings(dir.create(species_code))
  workdir <- species_code
  
  # Further subset sub_meta again
  sub_meta_species_code <- sub_meta[sub_meta$species_codes == species_code,]
  species_inds <- sub_meta_species_code$Sample_name
  
  # # Filter for individuals and invariants here
  # message("Subsetting VCF for individuals and removing invariants")
  # inds_to_keep <- paste(sub_meta_species_code$"Sample_name",collapse = ",")
  # system(paste0("bin/bcftools view -s ",inds_to_keep," ",vcf_path," -c 1 -c 1:nonmajor -Oz > ",species_code,"/subset_inds.vcf.gz"))
  # species_vcf <- paste0(species_code,"/subset_inds.vcf.gz")
  
  # Linkage prune using our current individuals
  set.seed(1000)
  message("Pruning for linkage of r2 0.2 with SNPRelate")
  snpset <- snpgdsLDpruning(gds_input, sample.id = species_inds, ld.threshold=0.2,autosome.only = F,num.thread = n_cores,slide.max.bp=50000L)
  
  # Subset for our unlinked SNPs
  if(length(unlist(snpset)) > snp_downsample){
    random_unlinked <- sample(unlist(snpset),snp_downsample)
  } else {
    random_unlinked <- unlist(snpset)
  }
  
  # Perform PCA
  message("Performing PCA and Fst structure assessment")
  pca <- snpgdsPCA(gds_input, sample.id = species_inds, snp.id=random_unlinked, num.thread=n_cores,autosome.only = F,eigen.cnt = 0)
  
  # Primary eigenvector
  eig1 <- pca$varprop[1]
  
  # 50% variance
  eig_cumsum <- cumsum(pca$varprop)
  eig50 <- length(eig_cumsum[eig_cumsum < 0.5]) + 1
  
  # Save the PCA as an R object
  #saveRDS(pca,paste0(species_code,"/SNPRelate_pca_results.rds"))
  
  # Now do FST based on population assignments...
  popmap <- data.frame(ind=sub_meta_species_code$Sample_name,
                       pop=paste0(sub_meta_species_code$Lat,"_",sub_meta_species_code$Long))
  snprelate_inds <- read.gdsn(index.gdsn(gds_input, "sample.id"))
  snprelate_inds <- snprelate_inds[snprelate_inds %in% popmap$ind]
  popmap <- popmap[match(snprelate_inds,popmap$ind),]
  
  # Produce a figure of PCA
  pca_res <- data.frame(ind=pca$sample.id,
                        PC1=pca$eigenvect[,1],
                        PC2=pca$eigenvect[,2])
  pca_res$pop <- popmap$pop
  pca_fig <- ggplot(pca_res,aes(PC1,PC2,colour=pop))+
    geom_point()+
    theme_minimal()+
    theme(axis.text = element_text(size=14),
          axis.title = element_text(size=16))+
    labs(x=paste0("PC1 (",round(pca$varprop[1]*100,2),"%)"),
         y=paste0("PC2 (",round(pca$varprop[2]*100,2),"%)"))+
    ggtitle(gsub("outputs/","",species_code))
  
  # Save to results dir
  pdf(paste0(species_code,"/SNPRelate_pca_fig.pdf"),width=10,height=8)
  print(pca_fig)
  dev.off()
  
  # Caculate W+C Fst
  species_fst <- snpgdsFst(gds_input,sample.id = popmap$ind,snp.id = random_unlinked, population=as.factor(popmap$pop),method="W&C84",autosome.only = F)
  
  # Save everything as a list and an RDS
  saveRDS(list(pca=pca,fst=species_fst,eig1=eig1,eig50=eig50),paste0(species_code,"/SNPRelate_pca_fst_results.rds"))
  
  # Finally, estimate population numbers based on DAPC
  # genotypes <- snpgdsGetGeno(gds_input, snp.id = random_unlinked)
  # genind_obj <- new.gen
  # grp <- find.clusters(x, max.n.clust=40)
}

##########################################################################################################################################################################################
# system(paste0("bin/bcftools annotate --set-id +'%CHROM\\:%POS' ",species_code,"/subset_inds.vcf.gz -Oz > ",species_code,"/subset_inds_annot.vcf.gz"))
# system(paste0("rm -f ",species_code,"/subset_inds.vcf.gz"))
# 
# # Linkage prune
# message("Linkage pruning for ",species_code)
# system(paste0("bin/plink2 --vcf ", species_code,"/subset_inds_annot.vcf.gz --out ",species_code,"/species_LD_plink_out_pruned --indep-pairwise 50 5 0.2 --allow-extra-chr --threads ",n_cores))
# 
# # Final filtering
# system(paste0("bin/bcftools filter -i 'ID=@",species_code,"/species_LD_plink_out_pruned.prune.in' ", species_code,"/subset_inds_annot.vcf.gz -Oz  > ",species_code,"/species_LD_pruned.vcf.gz"))
# 
# # Tidy up dir
# system(paste0("rm -f ",species_code,"/subset_inds_annot.vcf.gz ",species_code,"/species_LD_plink_out_pruned.prune.in ",species_code,"/species_LD_plink_out_pruned.prune.out"))
# 
# # Read in VCF
# vcf <- read.vcfR(paste0(species_code,"/species_LD_pruned.vcf.gz"))
# 
# # Subset for downsample
# vcf2 <- vcf[sample(1:nrow(vcf@gt),snp_downsample),]
# 
# # Perform PCA
# pca_gen <- df2genind(geno_matrix2,ncode = 1,ploidy = 2,sep="/")
# x.gen<-tab(pca_gen,freq=TRUE,NA.method="mean")





