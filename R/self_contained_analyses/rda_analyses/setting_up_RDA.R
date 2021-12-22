## RDA to assign variance to climate, geography, and population structure
lib <- c("pbmcapply","vegan","data.table","ggplot2","Rfast","dplyr")
sapply(lib,library,character.only=T)

# Set variables
n_cores <- 4

# For now, fetch the kubota AFs
res_dir <- "outputs/GEA_res/Arabidopsis_halleri_Kubota_Individual"
sample_metadata <- read.csv("metadata/sample_species_vcf_author_map_v2_210817c.csv")

# Input is AFs, col per population
AFs <- readRDS(paste0(res_dir,"/pop_allele_frqs.rds"))
AF_mat <- matrix(nrow = nrow(AFs[[1]]),ncol=length(AFs))
for(i in 1:ncol(AF_mat)){
  AF_mat[,i] <- AFs[[i]]$ref_freq
}
rownames(AF_mat) <- paste0(AFs[[1]]$chr,"-",AFs[[1]]$bp)
colnames(AF_mat) <- sapply(AFs,function(x) unique(x$pop))

# Impute missing genotypes as the median...
# row_medians <- rowMedians(AF_mat,na.rm=T,parallel=FALSE)
impute_medians <- function(genotype_row){
  tmp_row <- genotype_row
  tmp_row[is.na(tmp_row)] <- median(tmp_row,na.rm=T)
  return(tmp_row)
}
AF_mat_impute <- apply(AF_mat,1,impute_medians)

# Prepare Climate Inputs --------------------------------------------------
climate_cline <- read.table(paste0(res_dir,"/climate_cline.tsv"),header=T)
climate_vars <- colnames(climate_cline)[!(colnames(climate_cline) %in% c("Lat","Long"))]
climate_cline$pop <- paste0(climate_cline$Lat,"_",climate_cline$Long)

# Remove any rows with missing climate data, nothing can be done for them...
pops_to_keep <- na.omit(climate_cline)$pop
AF_mat_trans <- AF_mat_impute[pops_to_keep,]

# Standardization of the variables
scaled_climate <- scale(climate_cline[,climate_vars], center=TRUE, scale=TRUE) # center=TRUE, scale=TRUE are the defaults for scale()

# Recovering scaling coefficients
scale_env <- attr(scaled_climate, 'scaled:scale')
center_env <- attr(scaled_climate, 'scaled:center')

# Climatic table
scaled_climate <- as.data.frame(scaled_climate)
row.names(scaled_climate) <- climate_cline$pop

pop_structure <- readRDS(paste0(res_dir,"/SNPRelate_pca_fst_results.rds"))
struct_PCA <- data.frame(pop_structure$pca$eigenvect[,1:3],
                         ind=pop_structure$pca$sample.id)
colnames(struct_PCA) <- c(paste0("PC",1:3),"ind")

# Prepare Pop Structure Inputs --------------------------------------------
# Add in the population ID from the sample metadata...
sub_meta <- sample_metadata[sample_metadata$Sample_name %in% struct_PCA$ind,]
sub_meta$pop <- paste0(sub_meta$Lat,"_",sub_meta$Long)
sub_meta$ind <- sub_meta$Sample_name
struct_PCA <- merge(struct_PCA,sub_meta[,c("ind","pop")],by="ind")

# Fetch pop avgs...
pop_avg_PCA <- data.frame(struct_PCA %>% group_by(pop) %>% summarise(PC1_mean=mean(PC1,na.rm=T),
                                                                     PC2_mean=mean(PC2,na.rm=T),
                                                                     PC3_mean=mean(PC3,na.rm=T)))
rownames(pop_avg_PCA) <- pop_avg_PCA$pop
pop_avg_PCA <- pop_avg_PCA[rownames(scaled_climate),]


# Merge all the inputs together -------------------------------------------
explanatory_merge <- data.frame(climate_cline[,c("Lat","Long")],pop_avg_PCA,scaled_climate)

# Variable selection with forward model building --------------------------
## Null model
RDA0 <- rda(AF_mat_trans ~ 1,  explanatory_merge) 

## Full model
RDAfull <- rda(AF_mat_trans ~ mean_temp + mean_diurnal + isothermality + temp_seasonality + max_temp_warmest_month + min_temp_coldest_month + temp_range + mean_temp_wet_quarter + mean_temp_dry_quarter + mean_temp_warm_quarter + mean_temp_cold_quarter + annual_precip + precip_wet_month + precip_dry_month + precip_seasonality + precip_wet_quarter + precip_dry_quarter + precip_warm_quarter + precip_cold_quarter,explanatory_merge)

## Stepwise procedure with ordiR2step function
mod <- ordiR2step(RDA0, RDAfull, Pin = 0.05, R2permutations = 1000, R2scope = T)
climate_vars_to_keep <- colnames(attr(mod$terms,"factors"))

# Plot a correlation of climate vars
library(pheatmap)
pheatmap(cor(explanatory_merge[,climate_vars]))

# Variance partitioning ---------------------------------------------------
## Full model
# pRDAfull <- rda(AF_mat_trans ~ PC1_mean + PC2_mean + PC3_mean + Long + Lat + ,explanatory_merge)
pRDAfull <- rda(AF_mat_trans,explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean","Long","Lat",climate_vars_to_keep)])
RsquareAdj(pRDAfull)
# anova(pRDAfull) # Can't run this without killing R

## Pure climate model
pRDAclim <- rda(X=AF_mat_trans,Y=explanatory_merge[,c(climate_vars_to_keep)],Z=explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean","Long","Lat")]) 

## Pure neutral population structure model  
pRDAstruct <- rda(X=AF_mat_trans,Y=explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean")],Z=explanatory_merge[,c("Long","Lat",climate_vars_to_keep)])
RsquareAdj(pRDAstruct)

## Pure geography model 
pRDAgeog <- rda(X=AF_mat_trans,Y=explanatory_merge[,c("Long","Lat")],Z=explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean",climate_vars_to_keep)])
RsquareAdj(pRDAgeog)

## Assemble summary table of variance explained by each...
variance_table <- matrix(ncol=4,nrow=7)
rownames(variance_table) <- c("Full Model","Climate Model","Structure Model","Geography Model","Confounded Climate/Structure/Geography","Total unexplained","Total inertia")
colnames(variance_table) <- c("Inertia","R-squared","Proportion of explainable variance","Proportion of total variance")

# Fill inertia
variance_table[1,1] <- round(pRDAfull$CCA$tot.chi,1)
variance_table[2,1] <- round(pRDAclim$CCA$tot.chi,1)
variance_table[3,1] <- round(pRDAstruct$CCA$tot.chi,1)
variance_table[4,1] <- round(pRDAgeog$CCA$tot.chi,1)
variance_table[5,1] <- variance_table[1,1] - sum(variance_table[2:4,1])
variance_table[6,1] <- round(pRDAfull$CA$tot.chi,1)
variance_table[7,1] <- round(pRDAfull$tot.chi,1)

# Fill R2
variance_table[1,2] <- round(unlist(RsquareAdj(pRDAfull)[2]),3)
variance_table[2,2] <- round(unlist(RsquareAdj(pRDAclim)[2]),3)
variance_table[3,2] <- round(unlist(RsquareAdj(pRDAstruct)[2]),3)
variance_table[4,2] <- round(unlist(RsquareAdj(pRDAgeog)[2]),3)

# Fill explainable variance...
variance_table[1,3] <- 1
variance_table[2,3] <- round(variance_table[2,1]/variance_table[1,1],3)
variance_table[3,3] <- round(variance_table[3,1]/variance_table[1,1],3)
variance_table[4,3] <- round(variance_table[4,1]/variance_table[1,1],3)
variance_table[5,3] <- round(variance_table[5,1]/variance_table[1,1],3)

# Fill proportion of total variance
variance_table[,4] <- round(variance_table[,1] / variance_table[7,1],3)

# Variance partitioning of individual climate variables ------------------------------------------
climate_var_confound <- pbmclapply(climate_vars,function(climate_var){
  
  # Build a "full model" of both PCA structure, geography and the singular climate variable
  # pRDAfull_tmp <- rda(X=AF_mat_trans,Y=explanatory_merge[,c(climate_var,"PC1_mean","PC2_mean","PC3_mean","Long","Lat")]) 
  
  # Build a "partial model" conditioning on population structure and geography
  pRDAclim_tmp <- rda(X=AF_mat_trans,Y=explanatory_merge[,c(climate_var)],Z=explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean","Long","Lat")]) 
  
  # Build a "partial model" conditioning only on population structure
  pRDAclim_tmp2 <- rda(X=AF_mat_trans,Y=explanatory_merge[,c(climate_var)],Z=explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean")]) 
  
  # List outputs
  return(list(climate_var=climate_var,
              total_variance=pRDAclim_tmp$tot.chi,
              geog_struct_confound=pRDAclim_tmp$pCCA$tot.chi/pRDAclim_tmp$tot.chi,
              geog_struct_climate=pRDAclim_tmp$CCA$tot.chi/pRDAclim_tmp$tot.chi,
              struct_confound=pRDAclim_tmp2$pCCA$tot.chi/pRDAclim_tmp$tot.chi,
              struct_climate=pRDAclim_tmp2$CCA$tot.chi/pRDAclim_tmp$tot.chi))
},mc.cores=n_cores)

# All variables variance matrix
climate_var_output <- matrix(nrow=4,ncol=length(climate_vars))
for(i in 1:ncol(climate_var_output)){
  climate_var_output[1,i] <- climate_var_confound[[i]]$geog_struct_confound
  climate_var_output[2,i] <- climate_var_confound[[i]]$geog_struct_climate
  climate_var_output[3,i] <- climate_var_confound[[i]]$struct_confound
  climate_var_output[4,i] <- climate_var_confound[[i]]$struct_climate
}
rownames(climate_var_output) <- c("Confounded by Structure + Geography","Pure Climate Expl.",
                                  "Confounded by Structure","Climate + Geography Expl.")
colnames(climate_var_output) <- climate_vars

# Visualise...
climate_confound_long <- reshape2::melt(climate_var_output)
climate_confound_long[climate_confound_long$Var1 %in% c("Pure Climate Expl.","Climate + Geography Expl."),] %>%
  ggplot(aes(x=Var2,y=value,fill=Var1))+
  geom_bar(stat="identity",position="dodge")+
  theme(axis.text.x = element_text(angle=45,hjust=1))
  
  
