### IBE vs IBD
lib <- c("ggplot2","data.table","dplyr")
sapply(lib,library,character.only=T)

# Functions...
calc_qcd <- function(input_vector,na.rm=T){
  if(na.rm){
    input_vector <- na.omit(input_vector)
  }
  quarts <- quantile(input_vector,probs=c(0.25,0.75))
  ((quarts[2]-quarts[1]))/sum(quarts)
}
calc_MAD_stan <- function(input_vector,na.rm=T){
  if(na.rm){
    input_vector <- na.omit(input_vector)
  }
  median(abs(median(input_vector)-input_vector))/median(input_vector)
}

# Read in all the RDA results
res_dirs <- list.files("outputs/GEA_res")
res_dirs <- na.omit(sapply(res_dirs,function(dir) ifelse(file.exists(paste0("outputs/GEA_res/",dir,"/variance_partitioning_RDA_results.rds")),dir,NA)))
res_dirs <- grep("Cardamine",res_dirs,invert = T,value = T)
rda_res <- lapply(res_dirs,function(dir){
  tmp <- readRDS(paste0("outputs/GEA_res/",dir,"/variance_partitioning_RDA_results.rds"))
})
names(rda_res)

# Compare climate vs geog R2
climate_vs_geog_R2 <- data.frame(geogR2 = unlist(lapply(rda_res,function(x) return(x["Geography R-squared",]))),
                                 climR2 = unlist(lapply(rda_res,function(x) return(x["Climate variable R-squared",]))),
                                 climate_var = unlist(lapply(rda_res,function(x) return(colnames(x)))),
                                 dataset = rep(names(rda_res),each=ncol(rda_res[[1]])))
ggplot(climate_vs_geog_R2,aes(x=geogR2,y=climR2,colour=dataset))+
  geom_point(show.legend = F)+
  facet_wrap(~climate_var)

# Compare climate vs geog var explained...
GSEC_vs_SPEC <- data.frame(SPEC = unlist(lapply(rda_res,function(x) return(x["Climate by Geog R-squared",]))),
                           GSEC = unlist(lapply(rda_res,function(x) return(x["Climate variable R-squared",]))),
                           climate_var = unlist(lapply(rda_res,function(x) return(colnames(x)))),
                           dataset = rep(names(rda_res),each=ncol(rda_res[[1]])))
ggplot(GSEC_vs_SPEC,aes(x=SPEC,y=GSEC,colour=dataset))+
  geom_point(show.legend = F)+
  facet_wrap(~climate_var)

# Plot pdfs per climate variable
pdf("figs/rough_per_climate_IBA_by_geogEnv.pdf",width=10,height=10)
for(climate_var in unique(climate_vs_geogbyenv$climate_var)){
  print(ggplot(climate_vs_geogbyenv[climate_vs_geogbyenv$climate_var==climate_var,],aes(x=geog_env_R2,y=climR2))+
          geom_point()+
          geom_text(aes(label=dataset))+
          ggtitle(climate_var))
}
dev.off()

# Multivariate assessment of effect of selection across datasets
climate_r2_pca_dd <- matrix(nrow=length(rda_res),ncol=ncol(rda_res[[1]]))
rownames(climate_r2_pca_dd) <- names(rda_res)
colnames(climate_r2_pca_dd) <- colnames(rda_res[[1]])
for(i in 1:nrow(climate_r2_pca_dd)){
  climate_r2_pca_dd[i,] <- scale(climate_vs_geog_R2[climate_vs_geog_R2$dataset==names(rda_res)[i],"climR2"])
}

climate_r2_pca <- prcomp(climate_r2_pca_dd,scale. = F,center = T)
climate_r2_pca_scores <- data.frame(climate_r2_pca$x[,1:10])
colnames(climate_r2_pca_scores) <- paste0("PC",1:10)
climate_r2_pca_scores$dataset <- rownames(climate_r2_pca_scores)

# Visualise
R2_PCA <- ggplot(climate_r2_pca_scores,aes(x=PC1,y=PC2))+
  geom_point()+
  geom_text(aes(label=dataset))+
  labs(x=paste0("PC1 = ",round((climate_r2_pca$sdev[1]^2/sum(climate_r2_pca$sdev^2))*100),"%"),
       y=paste0("PC2 = ",round((climate_r2_pca$sdev[2]^2/sum(climate_r2_pca$sdev^2))*100),"%"))

# Loadings on main PCs...
climate_r2_pca_loadings <- data.frame(climate_r2_pca$rotation[,1:2])
colnames(climate_r2_pca_loadings) <- paste0("PC",1:2)
climate_r2_pca_loadings$climate_var <- rownames(climate_r2_pca_loadings)
climate_r2_pca_loadings <- reshape2::melt(climate_r2_pca_loadings)
#

# Single measure of variable importance
climate_r2_pca_loadings2 <- climate_r2_pca_loadings
for(i in 1:ncol(climate_r2_pca_loadings)){
  climate_r2_pca_loadings[,i] <- abs(climate_r2_pca_loadings[,i]*(climate_r2_pca_loadings$sdev[i]^2/sum(climate_r2_pca_loadings$sdev^2)))
}
sort(rowSums(fitness_pca_loadings2))

R2_PCA_loadings <- ggplot(climate_r2_pca_loadings,aes(y=climate_var,x=value))+
  geom_bar(stat="identity")+
  facet_wrap(~variable,ncol=2)

# Build dendrogram
climate_dendro <- hclust(dist(climate_r2_pca_scores[,paste0("PC",1:10)]))
library(ggdendro)
R2_dendro <- ggdendrogram(climate_dendro, rotate = TRUE, size = 2)


# Compare against climate cline PCA ---------------------------------------
climate_clines_res <- data.frame(rbindlist(lapply(res_dirs,function(dir){
  tmp <- read.table(paste0("outputs/GEA_res/",dir,"/climate_cline.tsv"),header = T)
  tmp$dataset=dir
  return(tmp[,c(-1,-2)])
})))

climate_PCA <- prcomp(climate_clines_res[,colnames(climate_clines_res) != "dataset"],scale=T,center = T)
climate_PCA_scores <- data.frame(climate_PCA$x[,1:10])
colnames(climate_PCA_scores) <- paste0("PC",1:10)
climate_PCA_scores$dataset=climate_clines_res$dataset

# Avg through datasets
climate_PCA_scores_avg <- data.frame(climate_PCA_scores %>% group_by(dataset) %>% summarise(mean_PC1=mean(PC1),
                                                                                            mean_PC2=mean(PC2),
                                                                                            mean_PC3=mean(PC3),
                                                                                            mean_PC4=mean(PC4),
                                                                                            mean_PC5=mean(PC5),
                                                                                            mean_PC6=mean(PC6),
                                                                                            mean_PC7=mean(PC7),
                                                                                            mean_PC8=mean(PC8),
                                                                                            mean_PC9=mean(PC9),
                                                                                            mean_PC10=mean(PC10)))


# Build dendrogram
rownames(climate_PCA_scores_avg) <- climate_PCA_scores_avg$dataset
climate_dendro_clines <- hclust(dist(climate_PCA_scores_avg[,paste0("mean_PC",1:10)]))
library(ggdendro)
cline_dendro <- ggdendrogram(climate_dendro_clines, rotate = TRUE, size = 2)

# Plot both dendros
library(patchwork)
pdf("figs/rough_climateR2_PCA_figs.pdf",width=20,height=8)
cowplot::plot_grid(
  R2_dendro,
  R2_PCA,
  R2_PCA_loadings,
  ncol=3,rel_widths=c(1,2,1))
R2_dendro+ggtitle("Clim R2 Dendro") | cline_dendro+ggtitle("Clim Cline PCA Dendro")
dev.off()


# Merge GSEC with Peak Skew and compare -----------------------------------
skew_res <- readRDS("outputs/wza_skewness.rds")
skew_res$climate_dataset <- paste0(skew_res$climate_var,"-",skew_res$dataset)
climate_vs_geogbyenv$climate_dataset <- paste0(climate_vs_geogbyenv$climate_var,"-",climate_vs_geogbyenv$dataset)

# Merge
skew_res_merge <- merge(skew_res,climate_vs_geogbyenv[,c("climate_dataset","climR2")],by="climate_dataset")

# Plot
ggplot(skew_res_merge,aes(x=climR2,y=skew))+
  geom_point()+
  geom_smooth(method="lm")


# GSEC-SPEC-QCD -----------------------------------------------------------
all_climate_clines <- data.frame(rbindlist(lapply(unique(GSEC_vs_SPEC$dataset),function(dataset){
  tmp <- read.table(paste0("outputs/GEA_res/",dataset,"/climate_cline.tsv"),header=T)
  
  # Transform temperatures to kelvin
  temp_vars <- grep("temp",colnames(tmp))
  for(var in temp_vars){
    tmp[,var] <- (tmp[,var]/10)+273.15
  }
  tmp$dataset=dataset
  # out <- data.frame(dataset=dataset,
  #                   climate_var=colnames(tmp)[!(colnames(tmp) %in% c("Long","Lat"))])
  # out$qcd <- apply(tmp[,out$climate_var],2,calc_qcd)
  # return(out)
  return(tmp[,!(colnames(tmp) %in% c("Lat","Long"))])
})))

# global_qcd <- apply(all_climate_clines[,colnames(all_climate_clines) != "dataset"],2,calc_qcd)

# Trans to long
all_climate_clines_long <- reshape2::melt(all_climate_clines)
# dataset_qcd <- data.frame(all_climate_clines_long %>% 
#                             group_by(dataset,variable) %>%
#                             summarise(qcd=calc_qcd(value)))
dataset_qcd <- as.data.table(all_climate_clines_long)[,.(qcd=calc_qcd(.SD$value,na.rm = T),
                                                         MAD=calc_MAD_stan(.SD$value,na.rm = T)),by=.(dataset,variable)]
global_mad <- sapply(unique(dataset_qcd$variable),function(x) max(dataset_qcd[dataset_qcd$variable==x,"MAD"]))
names(global_mad) <- unique(dataset_qcd$variable)
# dataset_qcds <- data.frame(reshape2::melt(all_climate_clines) %>% group_by(dataset,variable) %>%
#                              summarise(qcd=calc_qcd(value)))
# dataset_standard_var <- data.frame(all_climate_clines_long %>% group_by(dataset,variable) %>%
#                              summarise(var_climate=var(value)))

# Standardise
for(i in 1:length(global_mad)){
  dataset_qcd[dataset_qcd$variable == names(global_mad)[i],"relative_MAD"] <- (dataset_qcd[dataset_qcd$variable == names(global_mad)[i],"MAD"]/global_mad[i])
}

# for(i in 1:length(global_var)){
#   dataset_standard_var[dataset_standard_var$variable == names(global_var)[i],"NBR"] <- log2(dataset_standard_var[dataset_standard_var$variable == names(global_var)[i],"var_climate"]/global_var[i])
# }

pop_counts <- data.frame(table(all_climate_clines$dataset))
pop_counts$dataset=pop_counts$Var1

dataset_qcd$climate_dataset <- paste0(dataset_qcd$variable,"-",dataset_qcd$dataset)
GSEC_vs_SPEC$climate_dataset <- paste0(GSEC_vs_SPEC$climate_var,"-",GSEC_vs_SPEC$dataset)
GSEC_vs_SPEC_merge <- merge(GSEC_vs_SPEC,dataset_qcd[,c("climate_dataset","relative_MAD","MAD")],by="climate_dataset")
GSEC_vs_SPEC_merge <- merge(GSEC_vs_SPEC_merge,pop_counts,by="dataset")

ggplot(GSEC_vs_SPEC_merge,aes(Freq,MAD))+
  geom_point()

# plot each
library(cowplot)
plot_grid(
  ggplot(GSEC_vs_SPEC_merge,aes(MAD,GSEC))+geom_point(),
  ggplot(GSEC_vs_SPEC_merge,aes(MAD,SPEC))+geom_point(),
  ncol=1
)
library("scatterplot3d")
scatterplot3d(GSEC_vs_SPEC_merge[,c("GSEC","SPEC","relative_MAD")], pch = 16, type="h")
cor(GSEC_vs_SPEC_merge[,c("GSEC","SPEC","rel_qcd")])

# Fetch static predictions
GSEC_vs_SPEC_merge[GSEC_vs_SPEC_merge$SPEC < 0.1 &
                     GSEC_vs_SPEC_merge$SPEC > 0 & 
                     GSEC_vs_SPEC_merge$rel_qcd < 0.1,]

plot_wza_manhattan("outputs/GEA_res/Eucalyptus_sideroxylon_Murray_Individual","mean_temp")

# What kind of datasets do we have...
landscape_pca <- prcomp(GSEC_vs_SPEC_merge[,c("GSEC","SPEC","rel_qcd")],scale. = T,center = T)
landscape_pca_scores <- data.frame(landscape_pca$x)
summary(landscape_pca)
landscape_pca$rotation
ggplot(landscape_pca_scores,aes(x=PC1,y=PC3))+
  geom_point()


ggplot(all_climate_clines_long[all_climate_clines_long$dataset=="Eucalyptus_sideroxylon_Murray_Individual" & 
                                                         all_climate_clines_long$variable== "precip_dry_month",],aes(x=value))+
  geom_histogram()

hist(all_climate_clines_long[all_climate_clines_long$dataset=="Eucalyptus_sideroxylon_Murray_Individual" & 
                          all_climate_clines_long$variable== "precip_dry_month","value"])
hist(all_climate_clines_long[
                               all_climate_clines_long$variable== "precip_dry_month","value"])


# Using RDA,  estimate fitness with SPEC + Geog/Clim models ---------------

# For each dataset we want to build the linear expectation for residual variance explained by climate given the correlation between this and SPEC. Then take residuals
fitness_residuals_list <- lapply(rda_res,function(rda_tmp){
  
  # Take the variables of interest
  # tmp <- data.frame(SPEC=rda_tmp["Climate by Geog R-squared",],
  #                   GSEC=rda_tmp["Climate variable R-squared",],
  #                   unconfounded_clim=rda_tmp["Climate Expl. (+ Struct)",]/sum(rda_tmp["Confounded by Geography",],rda_tmp["Climate Expl. (+ Struct)",]))
  # 
  
  tmp <- data.frame(SPEC=rda_tmp["Climate by Geog R-squared",],
                    GSEC=rda_tmp["Climate variable R-squared",],
                    Y=rda_tmp["Climate Expl. (+ Struct)",],
                    YbyZ=rda_tmp["Confounded by Geography",])
  
  # Calculate proportion of GSEC that is Y
  tmp$Y_prop <- tmp$Y/rowSums(tmp[,c("Y","YbyZ")])
  
  
  ggplot(tmp,aes(x=SPEC,y=unconfounded_clim))+
    geom_point()+
    geom_smooth(method="lm")
  
  # # Firstly, just quantify correlation coefficient among these...
  # tmp_corr <- cor.test(tmp$SPEC,tmp$unconfounded_clim,method = "kendall")
  # tmp_corr2 <- cor.test(tmp$SPEC,tmp$GSEC,method = "kendall")
  # 
  # # Now build linear model
  # lm_tmp <- glm(unconfounded_clim~SPEC,tmp,family="gaussian")
  # tmp$predicted_unconfounded <- predict(lm_tmp)
  # 
  # # Fetch residuals...
  # tmp$residual_fitness <- tmp$unconfounded_clim-tmp$predicted_unconfounded
  # 
  # # Also standardise these within the dataset...
  # tmp$residual_fitness_std <- scale(tmp$residual_fitness)
  # tmp$climate_var <- rownames(tmp)
  return(list(dd=tmp))
  # return(list(dd=tmp,
  #             corr1=data.frame(tau=tmp_corr$estimate,pval=tmp_corr$p.value),
  #             corr2=data.frame(tau=tmp_corr2$estimate,pval=tmp_corr2$p.value)))
})

# Format corrs
data.frame(rbindlist(lapply(fitness_residuals_list,'[[',2)))

fitness_residuals_corrs$dataset <- names(rda_res)
ggplot(fitness_residuals_corrs,aes(x=dataset,y=tau))+
  geom_bar()

# Format
fitness_residuals <- data.frame(rbindlist(lapply(fitness_residuals_list,'[[',1)))
fitness_residuals$dataset <- rep(names(rda_res),each=nrow(fitness_residuals)/length(rda_res))

ggplot(fitness_residuals,aes(x=GSEC,YbyZ))+geom_point()

# What are all fitness residuals
hist(fitness_residuals$residual_fitness)
hist(fitness_residuals$residual_fitness_std)

# Build the PCA mentioned before to discuss similarity of selection among datasets...
residual_fitness_mat <- matrix(ncol=length(unique(fitness_residuals$climate_var)),nrow=length(unique(fitness_residuals$dataset)))
rownames(residual_fitness_mat) <- unique(fitness_residuals$dataset)
colnames(residual_fitness_mat) <- unique(fitness_residuals$climate_var)
for(i in 1:nrow(residual_fitness_mat)){
  for(j in 1:ncol(residual_fitness_mat)){
    residual_fitness_mat[i,j] <- fitness_residuals[fitness_residuals$dataset == rownames(residual_fitness_mat)[i] &
                                                     fitness_residuals$climate_var == colnames(residual_fitness_mat)[j],"residual_fitness_std"]
  }
}

# Make a PCA...
fitness_pca <- prcomp(residual_fitness_mat,scale. = F,center = T)
summary(fitness_pca)

fitness_pca_scores <- data.frame(fitness_pca$x)
fitness_pca_scores$dataset <- rownames(fitness_pca_scores)
# Visualise
R2_PCA <- ggplot(fitness_pca_scores,aes(x=PC1,y=PC2))+
  geom_point()+
  geom_text(aes(label=dataset))+
  labs(x=paste0("PC1 = ",round((fitness_pca$sdev[1]^2/sum(fitness_pca$sdev^2))*100),"%"),
       y=paste0("PC2 = ",round((fitness_pca$sdev[2]^2/sum(fitness_pca$sdev^2))*100),"%"))

# Loadings on main PCs...
fitness_pca_loadings <- data.frame(fitness_pca$rotation)
colnames(fitness_pca_loadings) <- paste0("PC",1:2)
fitness_pca_loadings$climate_var <- rownames(fitness_pca_loadings)
fitness_pca_loadings <- reshape2::melt(fitness_pca_loadings)

# Single measure of variable importance
fitness_pca_loadings2 <- fitness_pca_loadings
for(i in 1:ncol(fitness_pca_loadings2)){
  fitness_pca_loadings2[,i] <- abs(fitness_pca_loadings[,i]*(fitness_pca$sdev[i]^2/sum(fitness_pca$sdev^2)))
}
sort(rowSums(fitness_pca_loadings2))

fitness_pca_loadings_bar <- ggplot(fitness_pca_loadings,aes(y=climate_var,x=value))+
  geom_bar(stat="identity")+
  facet_wrap(~variable,ncol=2)

# Build dendrogram
fitness_dendro <- hclust(dist(fitness_pca_scores[,paste0("PC",1:10)]))
library(ggdendro)
fitness_dendro_fig <- ggdendrogram(fitness_dendro, rotate = TRUE, size = 2)

