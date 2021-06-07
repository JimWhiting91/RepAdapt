####################################################################################
# This script pulls in WZA res from outputs/ and compares with various dataset statistics
lib <- c("corrplot","EnvStats","parallel","data.table","tidyr","dplyr","vcfR","pbmcapply","ggplot2","ggridges")
dummy <- suppressPackageStartupMessages(lapply(lib,function(x){ 
  if (!require(x,character.only = T)){ 
    install.packages(x,repos='https://utstat.toronto.edu/cran/')
    library(x,character.only = T)
  } else {
    library(x,character.only = T)
  }
}))

# Fetch our metadata
metadata <- read.csv("metadata/sample_species_vcf_author_map_v2_210519.csv")

# Find all of our WZA results
dataset_dir <- list.files("outputs/")
########################################
# Remove these lines eventually
dataset_dir <- dataset_dir[dataset_dir != "snptable_convert"]
dataset_dir <- dataset_dir[dataset_dir != "Arabidopsis_lyrata_Willi_PoolSeq"]
########################################
wza_res <- lapply(dataset_dir,function(x) list.files(paste0("outputs/",x),pattern = "WZA"))
names(wza_res) <- dataset_dir

# Fetch the climate vars
climate_vars <- gsub("_WZA_pergene.tsv","",wza_res[[2]])

######################################################################################################
# Make a dataframe of "dataset" variables
popstructure_res <- data.frame(rbindlist(mclapply(dataset_dir,function(dataset){
  
  # Read in the res
  res_tmp <- readRDS(paste0("outputs/",dataset,"/SNPRelate_pca_fst_results.rds"))
  
  # Also count number of pops
  popN <- list.files(paste0("outputs/",dataset,"/"),pattern=".frq")
  
  # And SNP N
  snpN <- system(paste0("wc -l outputs/",dataset,"/",popN[[1]]),intern = T)
  snpN <- as.integer(strsplit(snpN," ")[[1]][1])
  
  # Set up out mat
  data.frame(dataset=dataset,
             indN=length(res_tmp$pca$sample.id),
             popN=length(popN),
             snpN=snpN,
             eig1=res_tmp$eig1,
             eig50=res_tmp$eig50)
  
},mc.cores=8)))

# Tidy up
popstructure_res$ind_per_pop <- popstructure_res$indN/popstructure_res$popN
popstructure_res$eig50_prop <- popstructure_res$eig50/popstructure_res$indN

######################################################################################################
##### Collect all the clines into a list #####
all_clines <- lapply(dataset_dir,function(dataset){
  # Read in the res
  cline_tmp <- read.table(paste0("outputs/",dataset,"/climate_cline.tsv"),header=T)
})
names(all_clines) <- dataset_dir

# Fetch the coefficient of variation
all_coef_var <- data.frame(rbindlist(lapply(dataset_dir,function(x){
  cv_tmp <- apply(na.omit(all_clines[[x]][,3:ncol(all_clines[[x]])]),2,cv)
  out <- data.frame(climate_var=names(cv_tmp),
                    cv=cv_tmp,
                    dataset=x)
})))

######################################################################################################
##### Associate clines and structure #####
structure_cline_corrs <- mclapply(dataset_dir,function(dataset){
  
  # Read in the res
  res_tmp <- readRDS(paste0("outputs/",dataset,"/SNPRelate_pca_fst_results.rds"))
  
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
  cline_tmp <- all_clines[[dataset]]
  cline_tmp$pop <- paste0(cline_tmp$Lat,"_",cline_tmp$Long)
  
  # Merge together
  clines_PC_merge <- na.omit(merge(pop_PCs,cline_tmp,"pop"))
  
  # Correlate
  corr_out <- cor(clines_PC_merge[,c("PC1","PC2")],clines_PC_merge[,climate_vars],method = "kendall")
  
},mc.cores=8)

# Plot all of the correlations to a pdf
pdf("figs/structure_vs_cline_correlations.pdf")
for(i in 1:length(structure_cline_corrs)){
  corrplot(structure_cline_corrs[[i]],title = dataset_dir[i])
}
dev.off()


######################################################################################################
# First, plot the distributions of all per-gene WZA scores across all datasets as ridges
all_ridges <- mclapply(climate_vars,function(var){
  
  # Fetch results
  dataset_res <- data.frame(rbindlist(lapply(dataset_dir,function(dataset){
    if(file.exists(paste0("outputs/",dataset,"/",var,"_WZA_pergene.tsv"))){
      tmp <- data.frame(fread(paste0("outputs/",dataset,"/",var,"_WZA_pergene.tsv")))
      tmp$dataset <- dataset
      return(tmp)
    } else {
      return(NULL)
    }
  })))
  
  # Plot distributions
  ridge_fig <- ggplot(dataset_res[dataset_res$weiZ != Inf,],aes(x=log(weiZ),y=dataset))+
    stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025,0.5,0.975), alpha = 0.7)+
    #stat_density_ridges(quantile_lines = TRUE, quantiles = 2)+
    theme_minimal()+
    labs(x="WZA (Weighted per-gene Z) [log]",y="Dataset")+
    ggtitle(var)
  return(ridge_fig)
},mc.cores=19)

# Print to a fig
pdf("figs/wza_density_ridges_across_datasets.pdf",height=12,width=8)
for(i in 1:length(all_ridges)){
  print(i)
  print(all_ridges[[i]])
}
dev.off()

######################################################################################################
# Repeat, but organise by variable instead
all_ridges <- mclapply(dataset_dir,function(dataset){
  
  # Fetch results
  dataset_res <- data.frame(rbindlist(lapply(wza_res[[dataset]],function(res){
    if(file.exists(paste0("outputs/",dataset,"/",res))){
      tmp <- data.frame(fread(paste0("outputs/",dataset,"/",res)))
      return(tmp)
    } else {
      return(NULL)
    }
  })))
  
  # Plot distributions
  ridge_fig <- ggplot(dataset_res[dataset_res$weiZ != Inf,],aes(x=log(weiZ),y=climate_var))+
    stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025,0.5,0.975), alpha = 0.7)+
    #stat_density_ridges(quantile_lines = TRUE, quantiles = 2)+
    theme_minimal()+
    labs(x="WZA (Weighted per-gene Z) [log]",y="Dataset")+
    ggtitle(dataset)
  return(ridge_fig)
},mc.cores=8)

# Print to a fig
pdf("figs/wza_density_ridges_across_variables.pdf",height=12,width=8)
for(i in 1:length(all_ridges)){
  print(i)
  print(all_ridges[[i]])
}
dev.off()

######################################################################################################
# Is variance among variables within datasets explained by any aspects of a dataset

# Calculate variance around median for all datasets
dataset_medians <- sapply(dataset_dir,function(dataset){
  
  # Fetch results
  dataset_res <- data.frame(rbindlist(lapply(wza_res[[dataset]],function(res){
    if(file.exists(paste0("outputs/",dataset,"/",res))){
      tmp <- data.frame(fread(paste0("outputs/",dataset,"/",res)))
      return(tmp)
    } else {
      return(NULL)
    }
  })))
  
  med_sum <- data.frame(dataset_res %>% group_by(climate_var) %>% summarise(med=median(weiZ)))
  return(var(med_sum$med))
})

# Any link with dataset features?
dataset_medians
cor(dataset_medians,popstructure_res[,c(2,3,4,5,7,8)],method = "spearman")
#cor.test(dataset_medians,popstructure_res[,c(2,3,4,5,7)],method = "spearman")

# Plot out the popN and eig1 effects
pdf("figs/median_variance_among_datasets_corrs.pdf")
plot(rank(dataset_medians)~popstructure_res$popN)
plot(rank(dataset_medians)~popstructure_res$eig1)
dev.off()



######################################################################################################
# For each weiZ distribution, examine variance plus median shift from 0
weiZ_summaries <- data.frame(rbindlist(lapply(dataset_dir,function(dataset){
  
  # Fetch results
  dataset_res <- data.frame(rbindlist(lapply(wza_res[[dataset]],function(res){
    if(file.exists(paste0("outputs/",dataset,"/",res))){
      tmp <- data.frame(fread(paste0("outputs/",dataset,"/",res)))
      return(tmp)
    } else {
      return(NULL)
    }
  })))
  
  # Variances
  climate_variances <- data.frame(dataset_res %>% group_by(climate_var) %>% summarise(weiZ_var=var(weiZ)))
  
  # Medians
  climate_median <- data.frame(dataset_res %>% group_by(climate_var) %>% summarise(weiZ_med=median(weiZ)))
  
  # Combine
  climate_variances$weiZ_median <- climate_median$weiZ_med
  climate_variances$dataset <- dataset
  return(climate_variances)
})))

# Merge and plot with cv
weiZ_summaries$group <- paste0(weiZ_summaries$dataset,":",weiZ_summaries$climate_var)
all_coef_var$group <- paste0(all_coef_var$dataset,":",all_coef_var$climate_var)
weiZ_summaries_merge <- merge(weiZ_summaries,all_coef_var[,c("cv","group")],"group")

# Also merge with PC1 and PC2 structure corrs
all_structure_corrs <- data.frame(rbindlist(lapply(1:length(structure_cline_corrs),function(i){
  
  # Transform
  corr_dd <- data.frame(t(structure_cline_corrs[[i]]))
  corr_dd$climate_var <- rownames(corr_dd)
  corr_dd$dataset <- dataset_dir[i]
  corr_dd$group <- paste0(corr_dd$dataset,":",corr_dd$climate_var)
  
  return(corr_dd[,c("PC1","PC2","group")])
})))

weiZ_summaries_merge <- merge(weiZ_summaries_merge,all_structure_corrs,"group")


# Visualise
pdf("figs/wza_distributions_vs_climate_variation.pdf")
ggplot(weiZ_summaries_merge,aes(log(cv),log(weiZ_median),colour=dataset))+
  geom_point()+
  theme_minimal()+
  labs(x="Coefficient of Variation",y="Median of WZA")
ggplot(weiZ_summaries_merge,aes(log(cv),log(weiZ_var),colour=dataset))+
  geom_point()+
  theme_minimal()+
  labs(x="Coefficient of Variation",y="Variance of WZA")
ggplot(weiZ_summaries_merge,aes(dataset,log(weiZ_var)))+
  #geom_point()+
  geom_boxplot()+
  theme_minimal()+
  theme(axis.text.x=element_text(angle=90,hjust=1))
  labs(x="Dataset",y="Variance of WZA")
dev.off()

# Visualise
pdf("figs/wza_distributions_vs_structure.pdf")
ggplot(weiZ_summaries_merge,aes(abs(PC1),log(weiZ_median),colour=dataset))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_minimal()+
  labs(x="Corr with structure PC1",y="Median of WZA")
ggplot(weiZ_summaries_merge,aes(abs(PC1),log(weiZ_var),colour=dataset))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_minimal()+
  labs(x="Corr with structure PC1",y="Variance of WZA")
dev.off()

# and correlations
cor.test(log(weiZ_summaries_merge$cv),log(weiZ_summaries_merge$weiZ_median))
cor.test(log(weiZ_summaries_merge$cv),log(weiZ_summaries_merge$weiZ_var))

######################################################################################################
#### Magnitude of structure vs variance effect by dataset #####
dataset_slopes <- data.frame(rbindlist(lapply(dataset_dir,function(dataset){
  
  # Fetch slope
  mod1 <- lm(log(weiZ_var)~abs(PC1),weiZ_summaries_merge[weiZ_summaries_merge$dataset==dataset,])
  
  return(data.frame(dataset=dataset,slope=mod1$coefficients[2]))
})))

# Merge
slope_merge <- merge(popstructure_res,dataset_slopes,"dataset")

# Plot
pdf("figs/effect_of_structure_variance_slopes.pdf")
ggplot(slope_merge,aes(eig1,slope))+
  geom_point()+
  geom_smooth(method="lm")
dev.off()


