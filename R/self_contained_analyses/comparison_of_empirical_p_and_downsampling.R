# Comparison of WZA distributions, empirical pvals and downsampling
lib <- c("data.table","ggplot2","ggridges","dplyr","corrplot","ggupset","UpSetR")
lapply(lib,library,character.only=T)

# Analysis looks at Amaranthus
all_res <- list.files("outputs/GEA_res/Amaranthus_tuberculatus_Wright_Individual/",pattern = "WZA")
climate_vars <- unique(sapply(strsplit(all_res,"_WZA_"),'[[',1))

# Data groups
data_groups <- c("WZA","WZA_empiricalp_full","WZA_empiricalp_half","WZA_empiricalp_tenth")

# First produce plots of distributions for all..
wza_res <- lapply(climate_vars,function(var){
  
  var_res <- data.frame(rbindlist(lapply(data_groups,function(group){
    tmp <- data.frame(fread(paste0("outputs/GEA_res/Amaranthus_tuberculatus_Wright_Individual/",var,"_",group,"_pergene.tsv")))
    tmp$data_group <- group
    tmp
  })))
var_res
})

# Produce ggridges for all climate vars...
data.frame(rbindlist(wza_res)) %>% 
ggplot(aes(x=weiZ,y=data_group))+
  facet_wrap(~climate_var)+
  geom_density_ridges()

# what's the median and variance of all distributions?
distribution_summaries <- data.frame(data.frame(rbindlist(wza_res)) %>% 
  group_by(climate_var,data_group) %>%
  summarise(median=median(weiZ),
            var=var(weiZ)))

ggplot(distribution_summaries,aes(x=data_group,y=median))+
  geom_boxplot()+
  facet_wrap(~climate_var)+
  theme(axis.text.x = element_text(angle=45,hjust=1))

ggplot(distribution_summaries,aes(x=data_group,y=var))+
  geom_boxplot()+
  facet_wrap(~climate_var)+
  theme(axis.text.x = element_text(angle=45,hjust=1))

#### Association between variance/median and pop_structure #####
pop_struc <- readRDS("outputs/GEA_res/Amaranthus_tuberculatus_Wright_Individual/SNPRelate_pca_fst_results.rds")
metadata  <- read.csv("metadata/sample_species_vcf_author_map_v2_210519.csv")

# Remake climate corrs
dataset_dir <- c("Amaranthus_tuberculatus_Wright_Individual")
structure_cline_corrs <- lapply(dataset_dir,function(dataset){
  
  # Read in the res
  res_tmp <- readRDS(paste0("outputs/GEA_res/",dataset,"/SNPRelate_pca_fst_results.rds"))
  
  # Organise individuals to their geopops
  inds <- res_tmp$pca$sample.id
  ind_pop <- sapply(inds,function(ind){
    paste0(na.omit(metadata[metadata$Sample_name == ind,"Lat"]),"_",na.omit(metadata[metadata$Sample_name == ind,"Long"]))
  })
  
  # Build PCA res with pop
  pca_scores <- data.frame(PC1=res_tmp$pca$eigenvect[,1],
                           PC2=res_tmp$pca$eigenvect[,2],
                           ind=inds,
                           pop=ind_pop)
  
  # Fetch centroids
  pop_PCs <- data.frame(pca_scores %>% group_by(pop) %>% summarise(PC1=mean(PC1),
                                                                   PC2=mean(PC2)))
  
  # Add to cline
  cline_tmp <- read.table(paste0("outputs/GEA_res/",dataset,"/climate_cline.tsv"),header=T)
  cline_tmp$pop <- paste0(cline_tmp$Lat,"_",cline_tmp$Long)
  
  # Merge together
  clines_PC_merge <- na.omit(merge(pop_PCs,cline_tmp,"pop"))
  
  # Correlate
  corr_out <- reshape2::melt(cor(clines_PC_merge[,c("PC1","PC2")],clines_PC_merge[,climate_vars],method = "kendall"))
  colnames(corr_out) <- c("PC","climate_var","corr")
  corr_out$dataset <- dataset
  
  return(corr_out)
})

# Combine the PC1 corrs with the distribution summaries
pc1_corr <- structure_cline_corrs[[1]][structure_cline_corrs[[1]]$PC=="PC1",]

pc1_corr_merge <- merge(distribution_summaries,pc1_corr,by = "climate_var")
pc1_corr_merge <- pc1_corr_merge[pc1_corr_merge$data_group %in% c("WZA","WZA_empiricalp_full"),]

# Visualise
ggplot(pc1_corr_merge,aes(x=corr,y=median))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~data_group)
ggplot(pc1_corr_merge,aes(x=corr,y=var))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~data_group)

################################################################

# What's the correlation between per gene scores among groups?
group_corrs <- lapply(wza_res,function(x){
  
  # genes to keep
  counts <- table(x$gene_id)
  genes_to_keep <- names(counts[counts==4])
  wza_sub <- x[x$gene_id %in% genes_to_keep,]
  
  # Make cormat
  to_corr <- matrix(ncol=4,nrow=length(genes_to_keep))
  for(i in 1:ncol(to_corr)){
    to_corr[,i] <- wza_sub[wza_sub$data_group==data_groups[i],"weiZ"]
  }
  colnames(to_corr) <- data_groups
  corr_res <- cor(to_corr)
  corrplot(corr_res,method = "ellipse")
})

# Pull out the top genes from all and assess the overlap
top_gene_overlaps <- lapply(wza_res,function(x){
  
  # For each, pull the top 5%
  top5 <- lapply(data_groups,function(y){
    x_sub <- x[x$data_group == y,]
    x_sub[x_sub$weiZ > quantile(x_sub$weiZ,0.05),]
  })
  names(top5) <- data_groups
  
  # Make upset
  data.frame(rbindlist(top5)) %>%
    group_by(gene_id) %>%
    summarise(data_group=list(data_group)) %>%
    ggplot(aes(x=data_group)) +
    geom_bar() +
    scale_x_upset()
})

# Plot these all
cowplot::plot_grid(plotlist = top_gene_overlaps[1:5])
