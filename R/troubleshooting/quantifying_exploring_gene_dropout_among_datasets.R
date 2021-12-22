# Pairwise orthology based on OF2 Blast results...
lib <- c("qvalue","wrswoR","ggExtra","pbmcapply","parallel","data.table","ggplot2","viridis","ggridges","dplyr","patchwork","readr")
sapply(lib,library,character.only=T)

# Where are the blast results?
blast_dir <- "outputs/orthology/Results_210913_17_genomes_with_conifers_noAA_filter/WorkingDirectory/"
blast_files <- list.files(blast_dir,pattern="Blast")

# Where is all the metadata...
vcf_genome_map <- read.table("metadata/vcf_genome_gff_210922_map.txt")
dataset_meta <- vcf_genome_map[,c("V8","V5")]
colnames(dataset_meta) <- c("dataset","genome")

# What are the OF codes for each of the genomes
genome_OF_codes <- data.frame(prot=c("Ahalleri.faa",
                                     "Alyrata.faa",
                                     "Athaliana.faa",
                                     "Atubercatus.faa",
                                     "Bstricta.faa",
                                     "Crubella.faa",
                                     "Egrandis.faa",
                                     "Hannuus.faa",
                                     "Mtruncatula.faa",
                                     "Pabies.faa",
                                     "Pdeltoides.faa",
                                     "Phallii.faa",
                                     "Pmenziesii.faa",
                                     "Ptaeda.faa",
                                     "Ptremula.faa",
                                     "Ptrichocarpa.faa",
                                     "Qpetraea.faa"))
genome_OF_codes$genome <- gsub(".faa","",genome_OF_codes$prot)
genome_OF_codes$OF_code <- (1:nrow(genome_OF_codes))-1

# Read in all of our liftovers from gff to OF codes
gff_OF_liftovers <- list.files("data/reference_genomes/",recursive = TRUE,pattern="_proteome_to_OF_id_map.txt")
gff_OF_liftovers <- lapply(genome_OF_codes$genome,function(x) grep(x,gff_OF_liftovers,value=T))
names(gff_OF_liftovers) <- genome_OF_codes$genome

# Let's converge... -------------------------------------------------------
focal_datasets <- list.files("outputs/GEA_res")
focal_climate <- "mean_temp"

# Remove any focal_datasets without full wza outputs...
focal_datasets <- na.omit(sapply(focal_datasets,function(x){
  if(file.exists(paste0("outputs/GEA_res/",x,"/",focal_climate,"_WZA_TC_allgenes.rds"))){
    return(x)
  } else {
    return(NA)
  }
}))

converge_matrix <- matrix(ncol=length(focal_datasets),nrow=length(focal_datasets))
rownames(converge_matrix) <- colnames(converge_matrix) <- focal_datasets

converge_corr_matrix <- converge_enrich_matrix <- converge_p_matrix <- converge_matrix

# To do this, we'll just go row by row and fill the matrix
Z_perms=100
wza_cutoffs=c(0.95,0.99)

# # Make our output matrices...
# converge_enrich_matrix_list <- lapply(1:length(wza_cutoffs),function(i) matrix(ncol=length(focal_datasets),nrow=length(focal_datasets)))
# names(converge_enrich_matrix_list) <- paste0("Cutoff = ",wza_cutoffs)

coverage_mat <- matrix(nrow=length(focal_datasets),ncol=2)
for(row in 1:nrow(converge_matrix)){
    # row_res <- pbmclapply(1:nrow(converge_matrix),function(row){
    
    # Fetch our focal dataset WZA results...
    dataset1 <- focal_datasets[row]
    genome1 <- na.omit(dataset_meta[dataset_meta$dataset == dataset1,"genome"])
    wza1 <- readRDS(paste0("outputs/GEA_res/",dataset1,"/",focal_climate,"_WZA_TC_allgenes.rds"))
    wza1$gea_gene <- wza1$gene_id
    
    # Identify consistent cutoffs
    wza1_cutoffs <- quantile(wza1$weiZ,probs=wza_cutoffs)
    
    # Merge these with the relevant OF_code...
    OF_codes1 <- read.table(paste0("data/reference_genomes/",gff_OF_liftovers[genome1]),header=T)
    if(colnames(OF_codes1)[3]=="seqid"){
      OF_codes1$gea_gene <- paste0(OF_codes1$seqid,":",OF_codes1$start,"-",OF_codes1$end)
    }
    OF_codes1$full_OF <- paste0(genome_OF_codes[genome_OF_codes$genome == genome1,"OF_code"],"_",OF_codes1$OF_ID)
    wza1_merge <- merge(wza1[,c("gea_gene","weiZ")],OF_codes1[,c("gea_gene","full_OF")])
    wza1_merge$qseqid <- wza1_merge$full_OF
    
print(paste0(dataset1,": ",round(nrow(wza1_merge)/nrow(OF_codes1),3)," proteome genes with wza"))
coverage_mat[row,1] <- dataset1
coverage_mat[row,2] <- round(nrow(wza1_merge)/nrow(OF_codes1),3)
}

# Save as a dataframe and re-order
coverage_mat <- data.frame(coverage_mat)
colnames(coverage_mat) <- c("dataset","prop_proteome_with_wza")
coverage_mat[order(coverage_mat$prop_proteome_with_wza),]

# Fetch the Weigel dataset to see snps
weigel_snps <- data.frame(fread("data/weigel_Athaliana_snps"))
colnames(weigel_snps) <- c("chr","pos")
weigel_snps$rand <- rnorm(nrow(weigel_snps),0,1)
weigel_snps$next_snp_dist <- c(weigel_snps$pos[2:nrow(weigel_snps)]-weigel_snps$pos[1:(nrow(weigel_snps)-1)],NA)

weigel_snps[weigel_snps$chr=="Athal_NC_003070.9" &
              weigel_snps$next_snp_dist > 0,] %>%
  ggplot(aes(y=next_snp_dist,x=pos))+
  geom_step()

quantile(weigel_snps[weigel_snps$chr=="Athal_NC_003070.9" &
              weigel_snps$next_snp_dist > 0,"next_snp_dist"],probs = 0.99) 

weigel_snps[weigel_snps$chr=="Athal_NC_003070.9" &
              weigel_snps$next_snp_dist > 0,] %>%
  ggplot(aes(y=rand,x=pos))+
  geom_point()+
  geom_segment(data=OF_codes1[OF_codes1$seqname=="Athal_NC_003070.9",],aes(y=10,yend=10,x=start,xend=end),colour="red2")

weigel_snps$pos2=weigel_snps$pos+1
weigel_snp_regions <- toGRanges(weigel_snps[weigel_snps$chr=="Athal_NC_003070.9",c("chr","pos","pos2")])
weigel_gene_regions <- toGRanges(OF_codes1[OF_codes1$seqname=="Athal_NC_003070.9",c("seqname","start","end")])
overlap <- regioneR::overlapRegions(weigel_gene_regions,weigel_snp_regions)
overlap$gene <- paste0(overlap$chr,":",overlap$startA,"-",overlap$endA)
length(unique(overlap$gene))

# Get SNP density
weigel_snps_test <- weigel_snps[weigel_snps$chr=="Athal_NC_003070.9",]
winds <- seq(0,max(weigel_snps_test$pos),100000)
winds2 <- winds+100000
density <- sapply(1:length(winds),function(x){
  return(nrow(weigel_snps_test[weigel_snps_test$pos >= winds[x] &
                                 weigel_snps_test$pos < winds2[x],]))
})
to_plot <- data.frame(wind_start = winds,
                      density=density)
hist(weigel_snps_test$pos)
ggplot(to_plot,aes(wind_start,density))+
  geom_point()+
  geom_segment(data=OF_codes1[OF_codes1$seqname=="Athal_NC_003070.9",],aes(y=5000,yend=5100,x=start,xend=end),colour="red2")

# Compare with roux snps
roux_snps <- data.frame(fread("data/roux_Athaliana_snps"))
hist(roux_snps$POS)
    