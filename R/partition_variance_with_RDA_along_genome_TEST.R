## RDA to assign variance to climate, geography, and population structure
lib <- c("R.utils","doParallel","pbmcapply","vegan","data.table","ggplot2","dplyr")
dummy <- suppressPackageStartupMessages(lapply(lib,function(x){ 
  if (!require(x,character.only = T)){ 
    install.packages(x,repos='https://utstat.toronto.edu/cran/')
    library(x,character.only = T)
  } else {
    library(x,character.only = T)
  }
}))

# Fetch command args ------------------------------------------------------
args <- commandArgs(asValues = T,excludeReserved = T)[-1]
print(args)

res_dir <- paste0("outputs/GEA_res/",args$dataset_dir)
metadata_path <- args$metadata # Needed for now, but remove once we have updated quantify_pop_structure to do pop avgs for us...
n_cores <- args$n_cores
model_selection_pval <- args$model_selection_pval
is.pool <- args$pool
snp_downsample <- args$snp_downsample

if(any(is.null(c(res_dir,metadata_path)))){
  stop("ERROR: No GEA results and/or METADATA provided")
}
if(is.null(n_cores)){
  n_cores <- 1
}
if(is.null(model_selection_pval)){
  model_selection_pval <- 0.05
}
if(is.null(is.pool)){
  is.pool <- FALSE
}
if(is.null(snp_downsample)){
  snp_downsample <- 50000
}

# Set a seed
set.seed(1000)

# Fetch metadata
sample_metadata <- read.csv(metadata_path)

############################################################################
# Set dummy variables
n_cores <- 4
is.pool <- TRUE
window_snps <- 200
############################################################################

# For now, fetch the kubota AFs
res_dir <- "outputs/GEA_res/Arabidopsis_thaliana_Gunther_PoolSeq//"
sample_metadata <- read.csv("metadata/sample_species_vcf_author_map_v2_210817c.csv")

# Set up parallel cores
cl <- makeCluster(n_cores, type="FORK")
registerDoParallel(cl)

# Kill cluster
on.exit(stopCluster(cl))

# Set seed for repeatability
set.seed(1000)

# Function lib ------------------------------------------------------------
impute_medians <- function(genotype_row){
  tmp_row <- genotype_row
  tmp_row[is.na(tmp_row)] <- median(tmp_row,na.rm=T)
  return(tmp_row)
}

# Prepare AFs -------------------------------------------------------------
message(paste0(">>> RDA modelling for ",args$dataset_dir,"
               "))
message(">>> Preparing Allele Frequency Matrix
        ")
# Input is AFs, col per population
AFs <- readRDS(paste0(res_dir,"/pop_allele_frqs.rds"))
AF_mat <- matrix(nrow = nrow(AFs[[1]]),ncol=length(AFs))
for(i in 1:ncol(AF_mat)){
  AF_mat[,i] <- AFs[[i]]$ref_freq
}
rownames(AF_mat) <- paste0(AFs[[1]]$chr,"-",AFs[[1]]$bp)
colnames(AF_mat) <- sapply(AFs,function(x) unique(x$pop))

# Impute missing genotypes as the median AND compartmentalise
message(paste0(">>> Building AF windows along the genome with ",window_snps," SNPs
        "))

# Fetch N SNPs per chromosome
chr_snp_counts <- table(AFs[[1]]$chr)
chr_to_keep <- names(chr_snp_counts)[chr_snp_counts > window_snps]
chr_winds <- lapply(chr_to_keep,function(chr) {
  winds <- seq(1,chr_snp_counts[chr],window_snps)
  winds2 <- winds+window_snps-1
  if(max(winds2)-max(winds) != window_snps){
    winds <- winds[1:(length(winds)-1)]
    winds2 <- winds+window_snps-1
  }
  return(winds)
})
names(chr_winds) <- chr_to_keep

# Build windows for all chroms
AF_mat_winds <- lapply(chr_to_keep,function(chr){
  message(paste0(">>> Starting windows for ",chr))
  window_names <- grep(paste0(chr,"-"),rownames(AF_mat),value=T)
  
  # Loop over windows and impute any missing rows...
  chr_AF_mat_winds <- foreach(i = 1:length(chr_winds[[chr]])) %dopar% {
    apply(AF_mat[window_names[chr_winds[[chr]][i]:(chr_winds[[chr]][i]+window_snps-1)],],1,impute_medians)
  }
  return(chr_AF_mat_winds)
})

# Remove windows that span chromosomes...




# Prepare Climate Inputs --------------------------------------------------
message(">>> Preparing climate data inputs
        ")
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

# Prepare Pop Structure Inputs --------------------------------------------
message(">>> Preparing population structure inputs (PCs 1-3)
        ")

# Slightly different formats depending on data type...
if(is.pool){
  pop_structure <- readRDS(paste0(res_dir,"/SNPRelate_pca_fst_results.rds"))
  struct_PCA <- data.frame(pop_structure$pca$x[,1:3],
                           ind=rownames(pop_structure$pca$x))
  colnames(struct_PCA) <- c(paste0("PC",1:3),"ind")
} else {
  pop_structure <- readRDS(paste0(res_dir,"/SNPRelate_pca_fst_results.rds"))
  struct_PCA <- data.frame(pop_structure$pca$eigenvect[,1:3],
                           ind=pop_structure$pca$sample.id)
  colnames(struct_PCA) <- c(paste0("PC",1:3),"ind")
}

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
# message(">>> Performing model selection on all climate variables (can take a while)...
#         ")
# RDA0 <- rda(AF_mat_trans ~ 1,  explanatory_merge)
# 
# ## Full model
# RDAfull <- rda(AF_mat_trans ~ mean_temp + mean_diurnal + isothermality + temp_seasonality + max_temp_warmest_month + min_temp_coldest_month + temp_range + mean_temp_wet_quarter + mean_temp_dry_quarter + mean_temp_warm_quarter + mean_temp_cold_quarter + annual_precip + precip_wet_month + precip_dry_month + precip_seasonality + precip_wet_quarter + precip_dry_quarter + precip_warm_quarter + precip_cold_quarter,explanatory_merge)
# 
# ## Stepwise procedure with ordiR2step function
# climate_mod_select <- ordiR2step(RDA0, RDAfull, Pin = model_selection_pval, R2permutations = 1000, R2scope = T)
# climate_vars_to_keep <- colnames(attr(climate_mod_select$terms,"factors"))
# 
# # Report
# if(!(is.null(climate_vars_to_keep))){
#   message(paste0(">>> For climate models, retaining: ",paste(climate_vars_to_keep,collapse = ",")))
# } else {
#   message(">>> For climate models, no variables retained... skipping general models...")
# }
# 
# # # Plot a correlation of climate vars
# # library(pheatmap)
# # pheatmap(cor(explanatory_merge[,climate_vars]))
# 
# # Variance partitioning ---------------------------------------------------
# ## Full model
# # pRDAfull <- rda(AF_mat_trans ~ PC1_mean + PC2_mean + PC3_mean + Long + Lat + ,explanatory_merge)
# if(!(is.null(climate_vars_to_keep))){
#   message(">>> Building full model...
#         ")
#   pRDAfull <- rda(AF_mat_trans,explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean","Long","Lat",climate_vars_to_keep)])
#   RsquareAdj(pRDAfull)
#   # anova(pRDAfull) # Can't run this without killing R
#   
#   ## Pure climate model
#   message(">>> Building climate model...
#         ")
#   pRDAclim <- rda(X=AF_mat_trans,Y=explanatory_merge[,c(climate_vars_to_keep)],Z=explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean","Long","Lat")]) 
#   
#   ## Pure neutral population structure model
#   message(">>> Building structure model...
#         ")
#   pRDAstruct <- rda(X=AF_mat_trans,Y=explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean")],Z=explanatory_merge[,c("Long","Lat",climate_vars_to_keep)])
#   RsquareAdj(pRDAstruct)
#   
#   ## Pure geography model
#   message(">>> Building geography model...
#         ")
#   pRDAgeog <- rda(X=AF_mat_trans,Y=explanatory_merge[,c("Long","Lat")],Z=explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean",climate_vars_to_keep)])
#   RsquareAdj(pRDAgeog)
#   
#   ## Assemble summary table of variance explained by each...
#   variance_table <- matrix(ncol=4,nrow=7)
#   rownames(variance_table) <- c("Full Model","Climate Model","Structure Model","Geography Model","Confounded Climate/Structure/Geography","Total unexplained","Total inertia")
#   colnames(variance_table) <- c("Inertia","R-squared","Proportion of explainable variance","Proportion of total variance")
#   
#   # Fill inertia
#   variance_table[1,1] <- round(pRDAfull$CCA$tot.chi,1)
#   variance_table[2,1] <- round(pRDAclim$CCA$tot.chi,1)
#   variance_table[3,1] <- round(pRDAstruct$CCA$tot.chi,1)
#   variance_table[4,1] <- round(pRDAgeog$CCA$tot.chi,1)
#   variance_table[5,1] <- variance_table[1,1] - sum(variance_table[2:4,1])
#   variance_table[6,1] <- round(pRDAfull$CA$tot.chi,1)
#   variance_table[7,1] <- round(pRDAfull$tot.chi,1)
#   
#   # Fill R2
#   variance_table[1,2] <- round(unlist(RsquareAdj(pRDAfull)[2]),3)
#   variance_table[2,2] <- round(unlist(RsquareAdj(pRDAclim)[2]),3)
#   variance_table[3,2] <- round(unlist(RsquareAdj(pRDAstruct)[2]),3)
#   variance_table[4,2] <- round(unlist(RsquareAdj(pRDAgeog)[2]),3)
#   
#   # Fill explainable variance...
#   variance_table[1,3] <- 1
#   variance_table[2,3] <- round(variance_table[2,1]/variance_table[1,1],3)
#   variance_table[3,3] <- round(variance_table[3,1]/variance_table[1,1],3)
#   variance_table[4,3] <- round(variance_table[4,1]/variance_table[1,1],3)
#   variance_table[5,3] <- round(variance_table[5,1]/variance_table[1,1],3)
#   
#   # Fill proportion of total variance
#   variance_table[,4] <- round(variance_table[,1] / variance_table[7,1],3)
# } else {
#   variance_table <- NULL
# }

# Variance partitioning of individual climate variables ------------------------------------------
# climate_var_confound <- pbmclapply(climate_vars,function(climate_var){
message(">>> Starting per-climate variable modelling
        ")
climate_var_confound <- foreach(climate_var = climate_vars) %dopar% {
  
  # First, just get the full R2 of this variable.
  RDA_var_tmp <- rda(X=AF_mat_trans,Y=explanatory_merge[,c(climate_var)])
  
  # Build a "partial model" conditioning on population structure and geography
  pRDAclim_tmp <- rda(X=AF_mat_trans,Y=explanatory_merge[,c(climate_var)],Z=explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean","Long","Lat")]) 
  
  # Build a "partial model" conditioning only on population structure
  pRDAclim_tmp2 <- rda(X=AF_mat_trans,Y=explanatory_merge[,c(climate_var)],Z=explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean")]) 
  
  # Build a "partial model" conditioning only on geog
  pRDAgeog_tmp <- rda(X=AF_mat_trans,Y=explanatory_merge[,c(climate_var)],Z=explanatory_merge[,c("Long","Lat")])
  
  # Reverse the model to condition geography on climate
  pRDAgeog_rev_tmp <- rda(X=AF_mat_trans,Z=explanatory_merge[,c(climate_var)],Y=explanatory_merge[,c("Long","Lat")])
  
  # Build a model of SNPs by geography
  GSSC_model_tmp <- rda(X=AF_mat_trans,Y=explanatory_merge[,c("Long","Lat")])
  
  # Build a model of climate var by geog...
  RDAgeog_clim_tmp <- rda(X=explanatory_merge[,c(climate_var)],Y=explanatory_merge[,c("Long","Lat")])
  
  # List outputs
  return(list(climate_var=climate_var,
              var_R2=RsquareAdj(RDA_var_tmp)$r.squared,
              var_R2_adj=RsquareAdj(RDA_var_tmp)$adj.r.squared,
              total_variance=pRDAclim_tmp$tot.chi,
              geog_struct_confound=pRDAclim_tmp$pCCA$tot.chi/pRDAclim_tmp$tot.chi,
              geog_struct_climate=pRDAclim_tmp$CCA$tot.chi/pRDAclim_tmp$tot.chi,
              struct_confound=pRDAclim_tmp2$pCCA$tot.chi/pRDAclim_tmp2$tot.chi,
              struct_climate=pRDAclim_tmp2$CCA$tot.chi/pRDAclim_tmp2$tot.chi,
              geog_confound=pRDAgeog_tmp$pCCA$tot.chi/pRDAgeog_tmp$tot.chi,
              geog_climate=pRDAgeog_tmp$CCA$tot.chi/pRDAgeog_tmp$tot.chi,
              geog_rev_confound=pRDAgeog_rev_tmp$pCCA$tot.chi/pRDAgeog_rev_tmp$tot.chi,
              geog_rev_climate=pRDAgeog_rev_tmp$CCA$tot.chi/pRDAgeog_rev_tmp$tot.chi,
              pure_climate=RDA_var_tmp$CCA$tot.chi/RDA_var_tmp$tot.chi,
              geog_climate_R2=RsquareAdj(RDAgeog_clim_tmp)$r.squared,
              geog_climate_R2_adj=RsquareAdj(RDAgeog_clim_tmp)$adj.r.squared)
  )
  # },mc.cores=n_cores)
}

# Finally get IBD effect...
RDA_geog <- rda(X=AF_mat_trans,Y=explanatory_merge[,c("Lat","Long")])
RDA_struct <- rda(X=AF_mat_trans,Y=explanatory_merge[,c("PC1_mean","PC2_mean","PC3_mean")])

# All variables variance matrix
climate_var_output <- matrix(nrow=13,ncol=length(climate_vars))
for(i in 1:ncol(climate_var_output)){
  climate_var_output[1,i] <- climate_var_confound[[i]]$geog_struct_confound
  climate_var_output[2,i] <- climate_var_confound[[i]]$geog_struct_climate
  climate_var_output[3,i] <- climate_var_confound[[i]]$struct_confound
  climate_var_output[4,i] <- climate_var_confound[[i]]$struct_climate
  climate_var_output[5,i] <- climate_var_confound[[i]]$geog_confound
  climate_var_output[6,i] <- climate_var_confound[[i]]$geog_climate
  climate_var_output[7,i] <- climate_var_confound[[i]]$var_R2
  climate_var_output[8,i] <- RsquareAdj(RDA_geog)$adj.r.squared
  climate_var_output[9,i] <- RsquareAdj(RDA_struct)$adj.r.squared
  climate_var_output[10,i] <- climate_var_confound[[i]]$pure_climate
  climate_var_output[11,i] <- climate_var_confound[[i]]$geog_climate_R2
  climate_var_output[12,i] <- climate_var_confound[[i]]$geog_rev_confound
  climate_var_output[13,i] <- climate_var_confound[[i]]$geog_rev_climate
}
rownames(climate_var_output) <- c("Confounded by Structure + Geography",
                                  "Climate Expl.",
                                  "Confounded by Structure",
                                  "Climate Expl. (+ Geog)",
                                  "Confounded by Geography",
                                  "Climate Expl. (+ Struct)",
                                  "Climate variable R-squared",
                                  "Geography R-squared",
                                  "Structure R-squared",
                                  "Simple Climate-Model Expl.",
                                  "Climate by Geog R-squared",
                                  "Geog Confounded by Climate",
                                  "Geog Expl. (- Climate)")
colnames(climate_var_output) <- climate_vars

# # Visualise...
# climate_confound_long <- reshape2::melt(climate_var_output)
# climate_confound_long[climate_confound_long$Var1 %in% c("Pure Climate Expl.","Climate + Geography Expl."),] %>%
#   ggplot(aes(x=Var2,y=value,fill=Var1))+
#   geom_bar(stat="identity",position="dodge")+
#   theme(axis.text.x = element_text(angle=45,hjust=1))

# Save everything to an object
message(">>> Saving results and finishing...")
# saveRDS(list(variable_selection=climate_mod_select,
#              full_variance_table=variance_table,
#              per_variable_table=climate_var_output),
#         paste0(res_dir,"/variance_partitioning_RDA_results.rds"))
saveRDS(climate_var_output,
        paste0(res_dir,"/variance_partitioning_RDA_results.rds"))


# Test Gunther ------------------------------------------------------------
test <- t(climate_var_output)
new_cols <- gsub(" ","_",colnames(test))
test <- data.frame(test)
colnames(test) <- new_cols
test$SPEC <- test$`Climate_by_Geog_R-squared`
test$GSEC <- test$`Simple_Climate-Model_Expl.`
test$var_a <- test$Confounded_by_Geography
test$var_b <- test$`Climate_Expl._(+_Struct)`
test$var_c <- test$Geog_Confounded_by_Climate
test$var_d <- test$`Geog_Expl._(-_Climate)`

test$b_rank <- rank(test$var_b)

# Calculate mys as a+b-b+d
test$myst <- test$var_a - test$var_d

# Calculate
test$myst_uncertainty <- test$myst/test$GSEC

ggplot(test,aes(SPEC,GSEC))+
  geom_point()
ggplot(test,aes(SPEC,myst_uncertainty))+
  geom_point()

# # Test assumptions of SPEC vs Mys -----------------------------------------
# test <- t(climate_var_output)
# new_cols <- gsub(" ","_",colnames(test))
# test <- data.frame(test)
# colnames(test) <- new_cols
# test$SPEC <- test$`Climate_by_Geog_R-squared`
# test$GSEC <- test$`Simple_Climate-Model_Expl.`
# test$var_a <- test$Confounded_by_Geography
# test$var_b <- test$`Climate_Expl._(+_Struct)`
# test$var_c <- test$Geog_Confounded_by_Climate
# test$var_d <- test$`Geog_Expl._(-_Climate)`
# 
# # Calculate mys as a+b-b+d
# test$mys <- test$var_a - test$var_d
# 
# # Calculate the suspected f coefficient and f(mys)
# test$f_coef <- test$var_b/rowSums(test[,c("var_b","var_d")])
# test$f_mys <- test$mys*test$f_coef
# 
# # Calculate missing central venn
# test$mys_prop <- test$mys / (test$var_a + test$var_c - test$mys)
# 
# # Plots of mys
# plot(test$mys_prop~test$SPEC)
# cor.test(test$mys,test$SPEC)
# 
# plot((test$mys/rowSums(test[,c("var_b","var_a")]))~test$SPEC)
# cor.test((test$mys/rowSums(test[,c("var_b","var_a")])),test$SPEC)
# 
# # Plots of f-mys
# plot(test$f_mys~test$SPEC)
# cor.test(test$f_mys,test$SPEC)
# 
# plot((test$f_mys/rowSums(test[,c("var_b","f_mys")]))~test$SPEC)
# cor.test((test$f_mys/rowSums(test[,c("var_b","f_mys")])),test$SPEC)
# 
# test$var=rownames(test)
# test$GSEC_confident <- test$var_b/rowSums(test[,c("var_b","var_c")])
# ggplot(test,aes(y=var,x=GSEC))+
#   geom_bar(stat="identity")
# 
# ggplot(test,aes(y=GSEC,x=var_b))+
#   geom_point()
# 
# test[order(-test$SPEC),]
# test$weirdness <- test$SPEC/test$mys
# test[order(-test$weirdness),]
# 

