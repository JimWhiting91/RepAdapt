### IBE vs IBD
lib <- c("ggplot2","data.table","dplyr","pbmcapply")
sapply(lib,library,character.only=T)


# Functions ---------------------------------------------------------------
# Quartile Coefficient of Dispersion
calc_qcd <- function(input_vector,na.rm=T){
  if(na.rm){
    input_vector <- na.omit(input_vector)
  }
  quarts <- quantile(input_vector,probs=c(0.25,0.75))
  ((quarts[2]-quarts[1]))/sum(quarts)
}

# Median Absolute Deviation, standardised by median
calc_MAD_stan <- function(input_vector,na.rm=T){
  if(na.rm){
    input_vector <- na.omit(input_vector)
  }
  median(abs(median(input_vector)-input_vector))/median(input_vector)
}

# PCA approach to summarise per-variable variance...



# -------------------------------------------------------------------------



# Read in all the RDA results
res_dirs <- list.files("outputs/GEA_res")
res_dirs <- na.omit(sapply(res_dirs,function(dir) ifelse(file.exists(paste0("outputs/GEA_res/",dir,"/variance_partitioning_RDA_results.rds")),dir,NA)))
res_dirs <- grep("Cardamine",res_dirs,invert = T,value = T)

rda_res <- data.frame(rbindlist(lapply(res_dirs,function(dir){
  # Read in 
  tmp <- readRDS(paste0("outputs/GEA_res/",dir,"/variance_partitioning_RDA_results.rds"))
  if(ncol(tmp) == 12){
    # tidy
    tmp$dataset <- dir
    return(tmp)
  } else {
    return(NULL)
  }
})))

# Make corrections
rda_res$SPEC_adj[rda_res$SPEC_adj < 0] <- 0
rda_res$GSEC_adj[rda_res$GSEC_adj < 0] <- 0

sample_sizes <- sapply(res_dirs,function(x) return(nrow(read.table(paste0("outputs/GEA_res/",x,"/climate_cline.tsv"),header=T))))
sample_sizes <- sort(sample_sizes)
rda_res$dataset_F <- factor(rda_res$dataset,levels=names(sample_sizes))

# Plot GSEC vs SPEC
ggplot(rda_res,aes(x=SPEC_adj,y=GSEC_adj,colour=dataset))+
  geom_point(show.legend = F)+
  geom_smooth(method="lm",show.legend = F)+
  # facet_wrap(~dataset,scales = "free")+
  theme_minimal()+
  labs(x="SPEC",y="GSEC")+
  xlim(0,1)+ylim(0,1)+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))
ggplot(rda_res,aes(x=SPEC_adj,y=GSEC_adj,colour=dataset))+
  geom_point(show.legend = F)+
  geom_smooth(method="lm",show.legend = F)+
  facet_wrap(~dataset,scales = "free")+
  theme_minimal()+
  labs(x="SPEC",y="GSEC")

# Plot corr per dataset
dataset_corrs <- as.data.table(rda_res)[,.(cor.test(.SD$GSEC_adj,.SD$SPEC_adj,method="kendall")$estimate),by=dataset]
ggplot(dataset_corrs,aes(y=dataset,x=V1))+
  geom_bar(stat="identity")+
  labs(x="SPEC-GSEC (K Tau)",y="Dataset")


# Plot myst
ggplot(rda_res,aes(x=SPEC_adj,y=clim_clean_abs,colour=dataset))+
  geom_point(show.legend = F)+
  geom_smooth(method="lm",show.legend = F)+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))+
  labs(x="SPEC",y="Confounded GSEC\n(Diff between a and b)")

# Plot SPEC vs uncertainty
ggplot(rda_res,aes(x=SPEC_adj,y=clim_clean_abs,colour=dataset))+
  geom_point(show.legend = F)+
  geom_smooth(method="lm",show.legend = F)+
  # ylim(0,1)+
  theme_minimal()+
  facet_wrap(~dataset,scales = "free")

# And corr coefficients...
dataset_corrs_myst <- as.data.table(rda_res)[,.(cor.test(.SD$clim_clean_abs,.SD$SPEC_adj,method="kendall")$estimate),by=dataset]
ggplot(dataset_corrs_myst,aes(y=dataset,x=V1))+
  geom_bar(stat="identity")+
  labs(x="SPEC-Uncertainty (K Tau)",y="Dataset")

# # Distribution of B
# hist(rda_res$clim_clean_abs)
# ggplot(rda_res,aes(y=climate_var,x=log(var_b)))+
#   geom_violin(draw_quantiles = T)
# 
# # Summed rank of b scores
# ggplot(rda_res,aes(y=climate_var,x=b_rank))+
#   geom_boxplot()
# b_rank_sums <- rda_res %>% group_by(climate_var) %>% summarise(b_rank_sum=sum(b_rank))
# 
# ggplot(rda_res,aes(y=myst,x=var_b))+
#   geom_point()+
#   geom_smooth(method="lm")+
#   facet_wrap(~dataset,scales = "free")
# 
# ggplot(rda_res,aes(y=myst,x=var_b,colour=dataset))+
#   geom_point(show.legend = F)+
#   geom_smooth(method="lm",show.legend = F)

# Plot all ranges of b-GSEC for each climate var
b_to_GSEC_figs <- lapply(unique(rda_res$climate_var),function(climate_var){
  
  tmp <- rda_res[rda_res$climate_var == climate_var,c("dataset","clim_clean","GSEC_adj")]
  
  ggplot(tmp,aes(y=dataset,x=clim_clean))+
    geom_point(size=3)+
    geom_point(aes(y=dataset,x=GSEC_adj),size=3,shape=17)+
    ggtitle(climate_var)+
    geom_segment(aes(y=dataset,yend=dataset,x=clim_clean,xend=GSEC_adj))+
    labs(x="% AF Variance Explained",y="Dataset")
  
})

ggplot(rda_res[rda_res$climate_var %in% c("mean_temp","annual_precip"),c("climate_var","dataset","clim_clean","GSEC_adj")],aes(y=dataset,x=clim_clean))+
  geom_point(size=3)+
  geom_point(aes(y=dataset,x=GSEC_adj),size=3,shape=17)+
  geom_segment(aes(y=dataset,yend=dataset,x=clim_clean,xend=GSEC_adj))+
  labs(x="% AF Variance Explained",y="Dataset")+
  facet_wrap(~climate_var,ncol=2)

# Can also visualise for a given dataset the relative importance of different variables in structuring genetic variation...
ggplot(rda_res[rda_res$dataset %in% grep("Helianthus",unique(rda_res$dataset),value=T),c("climate_var","dataset","clim_clean","GSEC_adj")],aes(y=climate_var,x=clim_clean))+
  geom_point(size=3)+
  geom_point(aes(y=climate_var,x=GSEC_adj),size=3,shape=17)+
  geom_segment(aes(y=climate_var,yend=climate_var,x=clim_clean,xend=GSEC_adj))+
  labs(x="% AF Variance Explained",y="Climate\nVariable")+
  facet_wrap(~dataset,nrow=1)


# Bring in measures of dispersion... --------------------------------------
# GSEC-SPEC-QCD -----------------------------------------------------------
all_climate_clines <- data.frame(rbindlist(lapply(unique(rda_res$dataset),function(dataset){
  tmp <- read.table(paste0("outputs/GEA_res/",dataset,"/climate_cline.tsv"),header=T)
  
  # Transform temperatures to kelvin
  temp_vars <- grep("temp",colnames(tmp))
  for(var in temp_vars){
    tmp[,var] <- (tmp[,var]/10)+273.15
  }
  tmp$dataset=dataset
  return(tmp[,!(colnames(tmp) %in% c("Lat","Long"))])
})))

# Trans to long
all_climate_clines_long <- reshape2::melt(all_climate_clines)
dataset_qcd <- as.data.table(all_climate_clines_long)[,.(qcd=calc_qcd(.SD$value,na.rm = T),
                                                         MAD=calc_MAD_stan(.SD$value,na.rm = T)),by=.(dataset,variable)]

# Standardise by max...
global_mad <- sapply(unique(dataset_qcd$variable),function(x) max(dataset_qcd[dataset_qcd$variable==x,"MAD"]))
names(global_mad) <- unique(dataset_qcd$variable)
for(i in 1:length(global_mad)){
  dataset_qcd[dataset_qcd$variable == names(global_mad)[i],"relative_MAD"] <- (dataset_qcd[dataset_qcd$variable == names(global_mad)[i],"MAD"]/global_mad[i])
}

# Standardise by global variance...
global_mad2 <- sapply(unique(dataset_qcd$variable),function(x) calc_MAD_stan(all_climate_clines_long[all_climate_clines_long$variable==x,"value"]))
names(global_mad2) <- unique(dataset_qcd$variable)
for(i in 1:length(global_mad2)){
  dataset_qcd[dataset_qcd$variable == names(global_mad2)[i],"relative_MAD2"] <- (dataset_qcd[dataset_qcd$variable == names(global_mad2)[i],"MAD"]/global_mad2[i])
}

# Merge with above estimates of GSEC, SPEC etc...
rda_res$id <- paste0(rda_res$dataset,"-",rda_res$climate_var)
dataset_qcd$id <- paste0(dataset_qcd$dataset,"-",dataset_qcd$variable)
rda_res_merge <- merge(rda_res,dataset_qcd[,c("id","relative_MAD","relative_MAD2")],by="id")

ggplot(rda_res_merge,aes(GSEC_adj,relative_MAD,colour=dataset))+
  geom_point(show.legend = F)+
  facet_wrap(~dataset,scales = "free")+
  geom_smooth(method="lm",show.legend = F)+
  labs(y="Relative Variance (MAD)",x="GSEC")+
  theme(axis.title=element_text(size=18))


# Case studies... ---------------------------------------------------------
# Within helianthus, 



# Using b+GSEC uncertainty, what is similarity ,  compare datasets ----------------------------------------------

# Draw a value for each climte-var for each dataset...
rda_res$b_mean <- rowMeans(rda_res[,c("clim_clean","GSEC_adj")])
rda_res$b_sd <- abs(rda_res$b_mean - rda_res$clim_clean)/2

# Now rearrange to pca wide-form
pca_mat2 <- matrix(ncol=length(unique(rda_res$climate_var)),nrow=length(unique(rda_res$dataset)))
rownames(pca_mat2) <- unique(rda_res$dataset)

# Scale row-wise so we're asking which are most important relatively among datasets...
for(i in 1:nrow(pca_mat2)){
  pca_mat2[unique(rda_res$dataset)[i],] <- scale(unlist(rda_res[rda_res$dataset == unique(rda_res$dataset)[i],"b_mean"]))
}
colnames(pca_mat2) <- unique(rda_res$climate_var)

# PCA
pca_fixed_tmp <- prcomp(pca_mat2,scale. = F,center = T)

# Loop over 100 times...
pca_loops <- data.frame(rbindlist(pbmclapply(1:100,function(iter){
  set.seed(iter)
  
  # Draw value
  env_importance <- as.data.table(rda_res)[,.(var_explained=rnorm(1,mean=.SD$b_mean,sd=.SD$b_sd)),by=.(dataset,climate_var)]
  
  pca_mat <- pca_mat2
  for(i in 1:nrow(pca_mat)){
    pca_mat[unique(env_importance$dataset)[i],] <- scale(unlist(env_importance[env_importance$dataset == unique(env_importance$dataset)[i],"var_explained"]))
  }
  
  # Project new variables into the space...
  pca_projected_tmp <- predict(object=pca_fixed_tmp,newdata = pca_mat)

  # Return scores
  score_out <- data.frame(pca_projected_tmp)
  
  # Also transform each
  score_out_trans <- score_out
  for(i in 1:ncol(score_out_trans)){
    score_out_trans[,i] <- score_out_trans[,i] * (pca_fixed_tmp$sdev^2/sum(pca_fixed_tmp$sdev^2))[i] 
  }
  colnames(score_out_trans) <- paste0("PC",1:ncol(score_out_trans),"_EigScaled")
  score_out <- cbind(score_out,score_out_trans)
  score_out$iter <- iter
  score_out$dataset <- rownames(score_out)
  

  return(score_out)
},mc.cores=4)))

# Plot...
centroids <- as.data.table(pca_loops)[,.(mean_PC1=mean(.SD$PC1),
                                         mean_PC2=mean(.SD$PC2),
                                         mean_PC3=mean(.SD$PC3),
                                         mean_PC4=mean(.SD$PC4),
                                         mean_PC5=mean(.SD$PC5),
                                         mean_PC6=mean(.SD$PC6),
                                         mean_PC7=mean(.SD$PC7),
                                         mean_PC8=mean(.SD$PC8),
                                         mean_PC9=mean(.SD$PC9),
                                         mean_PC10=mean(.SD$PC10)),by=.(dataset)]
ggplot(pca_loops,aes(x=PC1,y=PC2,colour=dataset))+
  geom_point(alpha=0.2)+
  geom_point(data=centroids,aes(mean_PC1,mean_PC2,colour=dataset),size=5)+
  stat_ellipse()

# Plot the distance between centroids...
rownames(centroids) <- centroids$dataset
centroids2 <- as.matrix(as.data.frame(centroids)[,grep("mean_",colnames(centroids),value=T)])
rownames(centroids2) <- centroids$dataset

centroid_dist <- as.matrix(dist(centroids2))
rownames(centroid_dist) <- rownames(centroids2) 
library(pheatmap)
pheatmap(centroid_dist)

# Also sum per variable the importance to get a single value...
pca_mat_ranks <- pca_mat2
for(i in 1:nrow(pca_mat_ranks)){
  pca_mat_ranks[i,] <- rank(pca_mat_ranks[i,])
}
climate_sums <- data.frame(climate_var = colnames(pca_mat_ranks),
                           climate_sum = colSums(pca_mat_ranks))
climate_sums$climate_var_F <- factor(climate_sums$climate_var,levels=climate_sums[order(climate_sums$climate_sum),"climate_var"])
ggplot(climate_sums,aes(y=climate_var_F,x=climate_sum))+
  geom_bar(stat="identity")+
  theme_minimal()+
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=16))+
  labs(y="",x="Summed Rank Importance")
                           


# Compare with skew res ---------------------------------------------------
skew_res <- readRDS("outputs/all_GEA_skew_results_tmp.rds")
skew_res$id <- paste0(skew_res$dataset,"-",skew_res$climate_var)
rda_res_merge <- merge(rda_res,skew_res[,c("id","skew")],by="id")

# Visualise vs uncertainty and GSEC
ggplot(rda_res_merge,aes(x=clim_confound,y=skew,colour=dataset))+
  geom_point(show.legend = F)+
  geom_smooth(method="lm")
  theme_minimal()+
  theme(axis.text = element_text(size=16),
        axis.title = element_text(size=18))+
  labs(x="Climate confound (b)",y="WZA Skew")
