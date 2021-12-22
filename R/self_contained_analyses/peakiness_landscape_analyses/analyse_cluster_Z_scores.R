# Analysis of how clustering quantification varies among datasets and among variables within datasets...
lib <- c("data.table","ggplot2","viridis","ggridges","dplyr")
sapply(lib,library,character.only=T)

# What metadata are we using?
metadata <- read.csv("metadata/sample_species_vcf_author_map_v2_210817c.csv")

# Function Library --------------------------------------------------------
plot_wza_manhattan <- function(wza_res_dir,climate_var){
  
  # Fetch the GEA res
  wza_res <- readRDS(paste0(wza_res_dir,"/",list.files(wza_res_dir,pattern=paste0(climate_var,"_WZA_TC_allgenes.rds"))))
  
  # Transform for plotting...
  wza_res$chr <- sapply(strsplit(wza_res$gene_id,":",),'[[',1)
  pos <- sapply(strsplit(wza_res$gene_id,":",),'[[',2)
  wza_res$start <- as.integer(sapply(strsplit(pos,"-",),'[[',1))
  wza_res$end <- as.integer(sapply(strsplit(pos,"-",),'[[',2))
  wza_res$mid <- rowMeans(wza_res[,c("start","end")])
  
  # Plot
  ggplot(wza_res,aes(mid,weiZ))+
    geom_point()+
    facet_grid(.~chr, scales = "free", space='free',switch="x")+
    theme()+
    theme(axis.text=element_text(size=16),
          axis.title = element_text(size=18),
          title = element_text(size=20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing.x=unit(0.1, "lines"),
          strip.text = element_text(angle=90,hjust=1,size=14))+
    geom_hline(yintercept = quantile(wza_res$weiZ,0.99),colour="red2")+
    labs(x="Chr/Position",y="WZA")
  
}

# Coefficient of variation
cv_calc <- function(data){
  sd(data,na.rm = T) / mean(data,na.rm = T) * 100
}

# Quartile Coefficient of Dispersion
qcd_calc <- function(data){
  quarts <- quantile(data, prob=c(.25,.75),na.rm = T)
  (quarts[2] - quarts[1])/sum(quarts)
}


# Prelim Plotting ---------------------------------------------------------

# Go and find all of the outputs and clustering results...
gea_res <- list.files("outputs/GEA_res/")
cluster_res <- na.omit(sapply(gea_res,function(x){
  if(file.exists(paste0("outputs/GEA_res/",x,"/clustering_Z_score_results.rds"))){
    return(paste0(x,"/clustering_Z_score_results.rds"))
  } else {
    return(NA)
  }
}))

# And remove gea_res that do not have cluster res
gea_res <- as.character(gsub("/clustering_Z_score_results.rds","",cluster_res))

# From these, now read in all of the Z scores and populate a matrix
# Just focus on wingspan-clustering for now...
all_cluster_res <- lapply(paste0("outputs/GEA_res/",cluster_res),readRDS)

# Build and populate matrix
cluster_mat <- matrix(nrow=nrow(all_cluster_res[[1]]),ncol=length(all_cluster_res))
rownames(cluster_mat) <- all_cluster_res[[1]]$climate_var
colnames(cluster_mat) <- gea_res
for(i in 1:ncol(cluster_mat)){
  for(var in rownames(cluster_mat)){
    if(length(all_cluster_res[[i]][all_cluster_res[[i]]$climate_var == var,"wingspan_obs_Z"]) != 0){
      cluster_mat[var,i] <- all_cluster_res[[i]][all_cluster_res[[i]]$climate_var == var,"wingspan_obs_Z"]
    }
  }
}

# Visualise...
long_cluster <- reshape2::melt(cluster_mat)
colnames(long_cluster) <- c("climate_var","dataset","cluster_Z")
ggplot(long_cluster,aes(climate_var,dataset,fill=cluster_Z))+
  geom_tile()+
  scale_fill_viridis(option="A")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1))


# # TEMP SPACE --------------------------------------------------------------
# pdf("figs/rough_all_mean_temp_manhattans.pdf",width=12,height=4)
# for(dir in res_dirs){
#   print(plot_wza_manhattan(paste0("outputs/GEA_res/",dir),"mean_temp"))
# }
# dev.off()
# 
# ###########################################################################

# Plot the most vs least clustered...
# Least
plot_wza_manhattan("outputs/GEA_res/Arabidopsis_thaliana_Roux_PoolSeq","annual_precip")
plot_wza_manhattan("outputs/GEA_res/Arabidopsis_thaliana_Weigel_WGS","mean_temp")

# Most - These seem to correspond with haploblocks identified in Todesco
plot_wza_manhattan("outputs/GEA_res/Helianthus_annuus_Todesco_Individual","mean_diurnal")
plot_wza_manhattan("outputs/GEA_res/Helianthus_annuus_Todesco_Individual","annual_precip")
plot_wza_manhattan("outputs/GEA_res/Helianthus_argophyllus_Todesco_Individual","temp_range")
plot_wza_manhattan("outputs/GEA_res/Helianthus_petiolaris_Todesco_Individual","mean_diurnal")
plot_wza_manhattan("outputs/GEA_res/Helianthus_petiolaris_Todesco_Individual","precip_wet_month")

# plot_wza_manhattan("outputs/GEA_res/Eucalyptus_albens_Murray_Individual/","mean_temp_wet_quarter")
plot_wza_manhattan("outputs/GEA_res/Arabidopsis_halleri_Kubota_Individual","annual_precip")

# Visualise datasets and climate vars as distributions 
# Distribution by dataset
ggplot(long_cluster,aes(y=dataset,x=cluster_Z))+
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2)+
  theme_minimal()

# Distribution by cliamte_var
ggplot(long_cluster,aes(y=climate_var,x=cluster_Z))+
  stat_density_ridges(quantile_lines = TRUE, quantiles = 2)+
  theme_minimal()

# How does clustering link with variability and pop structure -------------

# For each dataset, build info on correlation of major structure axis with climate var
all_climate_clines <- lapply(paste0("outputs/GEA_res/",gea_res,"/climate_cline.tsv"),read.table,header=T)
all_structure_res <- lapply(paste0("outputs/GEA_res/",gea_res,"/SNPRelate_pca_fst_results.rds"),readRDS)
names(all_climate_clines) <- gea_res
names(all_structure_res) <- gea_res


# Association with pop structure ------------------------------------------
# We first need to build up correlation matrices for each dataset of structure PC vs climate...
climate_pop_structure_corr <- data.frame(rbindlist(lapply(gea_res,function(dataset){
  
  print(dataset)
  
  cline_tmp <- all_climate_clines[[dataset]]
  pca_tmp <- all_structure_res[[dataset]]
  
  # Fetch PC1 and PC2 scores...
  if(length(pca_tmp) == 4){
    pc_scores <- data.frame(pca_tmp$pca$eigenvect[,1:2])
    colnames(pc_scores) <- c("PC1","PC2")
    pc_scores$sample <- pca_tmp$pca$sample.id
    eigenvals <- pca_tmp$pca$eigenval
    eigenvals <- eigenvals/sum(eigenvals)*100
  } else if (length(pca_tmp) == 5){
    pc_scores <- data.frame(pca_tmp$pca$x[,1:2])
    colnames(pc_scores) <- c("PC1","PC2")
    pc_scores$sample <- rownames(pca_tmp$pca$x)
    eigenvals <- pca_tmp$pca$sdev^2/sum(pca_tmp$pca$sdev^2)*100
  }
  
  # Get our pops from the metadata...
  sub_meta <- metadata[metadata$Sample_name %in% pc_scores$sample,]
  sub_meta$geo_pops <- paste0(sub_meta$Lat,"_",sub_meta$Long)
  sub_meta$sample <- sub_meta$Sample_name
  
  # Add these to pc scores
  pc_scores <- merge(pc_scores,sub_meta[,c("sample","geo_pops")],by="sample")
  
  # Group these
  pop_pcs <- data.frame(pc_scores %>% group_by(geo_pops) %>% summarise(mean_PC1=mean(PC1,na.rm=T),
                                                                       mean_PC2=mean(PC2,na.rm=T)))
  
  # Add pop IDs to cline...
  cline_tmp$geo_pops <- paste0(cline_tmp$Lat,"_",cline_tmp$Long)
  cline_tmp_merge <- na.omit(merge(cline_tmp,pop_pcs,by="geo_pops"))
  
  # And make corr mat
  climate_vars <- colnames(cline_tmp)
  climate_vars <- climate_vars[!(climate_vars %in% c("Lat","Long","geo_pops"))]
  climate_corrs <- cor(cline_tmp_merge[,climate_vars],cline_tmp_merge[,paste0("mean_PC",1:2)],method = "spearman")
  
  # Transform and export...
  out <- reshape2::melt(climate_corrs)
  colnames(out) <- c("climate_var","PC","spearman_rho")
  out$PC <- gsub("mean_","",out$PC)
  out$dataset <- dataset
  
  # Also add variance explained...
  out$str_explained <- NA
  for(i in 1:2){
    out[out$PC == paste0("PC",i),"str_explained"] <- eigenvals[i]
  }
  return(out)
})))

# Also fetch the eig50s
dataset_eig50 <- sort(sapply(gea_res,function(dataset){
  pca_tmp <- all_structure_res[[dataset]]
  
  # For eig50, need the proportion rather than absolute
  if(length(pca_tmp) == 4){
    eigenvector_N <- ncol(pca_tmp$pca$eigenvect)
  } else if(length(pca_tmp) == 5){
    eigenvector_N <- ncol(pca_tmp$pca$rotation)
  }
  return(pca_tmp$eig50/eigenvector_N)
}))

# Order datasets by those with strongest eig1...
eig1_orders <- unique(climate_pop_structure_corr[climate_pop_structure_corr$PC == "PC1",c("dataset","str_explained")])
for(dataset in eig1_orders$dataset){
  eig1_orders[eig1_orders$dataset == dataset,"eig50"] <- dataset_eig50[dataset]
}
eig1_orders <- eig1_orders[order(eig1_orders$eig50),]
long_cluster$dataset_F <- factor(long_cluster$dataset,levels = eig1_orders$dataset)

# Add these to long_cluster
long_cluster$climate_structure_corr <- NA
for(i in 1:nrow(long_cluster)){
  long_cluster$climate_structure_corr[i] <- climate_pop_structure_corr[climate_pop_structure_corr$dataset==long_cluster$dataset[i] &
                                                                         climate_pop_structure_corr$climate_var==long_cluster$climate_var[i] &
                                                                         climate_pop_structure_corr$PC == "PC1","spearman_rho"]
}
long_cluster <- merge(long_cluster,eig1_orders,by="dataset")

# visualise climate corrs...
ggplot(long_cluster,aes(climate_var,dataset,fill=abs(climate_structure_corr)))+
  geom_tile()+
  scale_fill_viridis(option="A")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1))

# Visualise climate corrs vs cluster Z
ggplot(long_cluster,aes(abs(climate_structure_corr),cluster_Z,colour=dataset))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~dataset_F)+
  theme(legend.position = "none")

# Visualise climate corrs vs cluster Z with scaling
ggplot(long_cluster,aes(abs(climate_structure_corr),cluster_Z,colour=dataset))+
  geom_point()+
  geom_smooth(method="lm")+
  #facet_wrap(~dataset_F)+
  theme(legend.position = "none")

# Visualise climate corrs by dataset
ggplot(long_cluster,aes(y=dataset_F,x=abs(climate_structure_corr),fill=str_explained))+
  geom_boxplot(outlier.colour = NA)+
  theme(legend.position = "none")+
  scale_fill_viridis(option="A")

# Visualise cluster_z by dataset
ggplot(long_cluster,aes(y=dataset_F,x=cluster_Z,fill=str_explained))+
  geom_boxplot(outlier.colour = NA)+
  theme(legend.position = "none")+
  scale_fill_viridis(option="A")




# Building models ---------------------------------------------------------
long_cluster$scaled_climate_pc1_corr <- abs(long_cluster$climate_structure_corr) * long_cluster$str_explained
glm1 <- glm(cluster_Z~dataset+climate_var+str_explained+abs(climate_structure_corr)+scaled_climate_pc1_corr,data = long_cluster)
step(glm1,test="F")

# Association with CV -----------------------------------------------------
# Add the coefficient of variation for each of these...
all_climate_clines_cv <- lapply(all_climate_clines,function(x){
  cv_mat <- matrix(ncol=1,nrow=19)
  rownames(cv_mat) <- colnames(x)[3:ncol(x)]
  for(i in 1:nrow(cv_mat)){
    #cv_mat[i,1] <- abs(cv_calc(x[,rownames(cv_mat)[i]]))
    cv_mat[i,1] <- qcd_calc(x[,rownames(cv_mat)[i]])
    #cv_mat[i,1] <- var(scale(na.omit(x[,rownames(cv_mat)[i]])))
  }
  return(cv_mat)
})

# Add these to long_cluster
long_cluster$climate_CV <- NA
for(i in 1:nrow(long_cluster)){
  long_cluster$climate_CV[i] <- all_climate_clines_cv[[long_cluster$dataset[i]]][long_cluster$climate_var[i],1]
}

# visualise how variable climates are within datasets...
ggplot(long_cluster,aes(climate_var,dataset,fill=log10(abs(climate_CV))))+
  geom_tile()+
  scale_fill_viridis(option="A")+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=45,hjust=1))

# Visualise
ggplot(long_cluster,aes(abs(climate_CV),cluster_Z,colour=dataset))+
  geom_point()+
  geom_smooth(method="lm")+
 # facet_wrap(~dataset,scales="free")+
  theme(legend.position = "none")
