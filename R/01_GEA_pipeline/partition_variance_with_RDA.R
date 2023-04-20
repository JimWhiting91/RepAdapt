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

vcf_path <- args$vcf
res_dir <- paste0("outputs/GEA_res/",args$dataset_dir)
n_cores <- as.integer(args$n_cores)
snp_downsample <- as.integer(args$snp_downsample)

message(paste0(">>> VCF = ",vcf_path))
message(paste0(">>> Results Dir = ",res_dir))
message(paste0(">>> Core N = ",n_cores))
message(paste0(">>> SNP Downsample N = ",snp_downsample))

if(any(is.null(c(vcf_path,res_dir)))){
  stop("ERROR: No VCF or Output directory provided...")
}
if(is.null(n_cores)){
  n_cores <- 1
}
if(is.null(snp_downsample)){
  snp_downsample <- 10000
}

# ############################################################################
# # Set dummy variables
# vcf_path <- "data/VCFs/19c_Emagnificata_Murray/murray_Emag_full_concatened.vcf.gz"
# n_cores <- 12
# snp_downsample <- 10000
# res_dir <- "outputs/GEA_res/230321_Eucalyptus_magnificata_Murray_Individual/"

# ############################################################################

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

# Let's only retain cases with missing data <10% and downsample from those
perSNP_missing = rowSums(is.na(AF_mat))
to_keep = which(perSNP_missing/ncol(AF_mat) < 0.1)
AF_mat_lowMiss = AF_mat[to_keep,]

# Impute missing genotypes as the median AND downsample
# Can we downsample without having to impute?
no_missing_rows <- which(perSNP_missing == 0)

# Preferentially shuffle SNPs with full data, then preferentially shuffle HQ SNPs, then settle for whatever HQ we have...
if(length(no_missing_rows) > snp_downsample){
  message(">>> Downsampling 'No Missing' SNPs (100% data) to ",snp_downsample," SNPs")
  AF_mat_impute <- t(AF_mat[sort(sample(no_missing_rows,snp_downsample)),])
} else if(nrow(AF_mat_lowMiss) > snp_downsample) {
  message(">>> Downsampling HQ SNPs (>90% data) to ",snp_downsample," SNPs")
  AF_mat_impute <- apply(AF_mat_lowMiss[sort(sample(1:nrow(AF_mat_lowMiss),snp_downsample)),],1,impute_medians)
} else {
  message(">>> Not enough data for downsample, retaining ",nrow(AF_mat_lowMiss)," SNPs")
  AF_mat_impute <- apply(AF_mat_lowMiss,1,impute_medians)
}

# Prepare Climate Inputs --------------------------------------------------
message(">>> Preparing climate data inputs
        ")
climate_cline <- read.table(paste0(res_dir,"/climate_cline.tsv"),header=T)
climate_cline$pop <- paste0(climate_cline$Lat,"_",climate_cline$Long)

# Remove the climate change variables...
climate_cline = climate_cline[,grep("clim_change",colnames(climate_cline),invert = T)]
climate_vars <- colnames(climate_cline)[!(colnames(climate_cline) %in% c("Lat","Long","pop"))]

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

# Merge these with pop identifiers
climate_cline$pop <- paste0(climate_cline$Lat,"_",climate_cline$Long)
explanatory_merge <- data.frame(climate_cline[,c("pop","Long","Lat")],scaled_climate)
colnames(explanatory_merge)[1] <- "pop"

# Make sure we prune any populations that don't have climate data...
to_keep <- explanatory_merge$pop[!(is.na(explanatory_merge$mean_temp))]

# Variance partitioning of individual climate variables ------------------------------------------
message(">>> Starting per-climate variable modelling
        ")

# Remove any climate variables where we just don't have any data, i.e. column of NAs
clim_remove = names(which(apply(explanatory_merge,2,function(x) sum(is.na(x))) == nrow(explanatory_merge)))
climate_vars = climate_vars[!climate_vars %in% clim_remove]

# Alternative variance partitioning workflow... ---------------------------
variance_partitioning_clean <- foreach(climate_var = climate_vars,.combine=rbind) %dopar% {
  
  # Variance partition with space
  part1 <- varpart(AF_mat_trans[to_keep,],as.matrix(explanatory_merge[to_keep,climate_var]),explanatory_merge[to_keep,c("Long","Lat")])
  
  # Use the standard RDA models to get our a,b,c partitions
  rda1 <- rda(AF_mat_trans[to_keep,],as.matrix(explanatory_merge[to_keep,climate_var]),explanatory_merge[to_keep,c("Long","Lat")])
  rda2 <- rda(AF_mat_trans[to_keep,],explanatory_merge[to_keep,c("Long","Lat")],as.matrix(explanatory_merge[to_keep,climate_var]))
  a = rda1$CCA$tot.chi / rda1$tot.chi
  bc = rda1$pCCA$tot.chi / rda1$tot.chi
  c = rda2$CCA$tot.chi / rda2$tot.chi
  ab = rda2$pCCA$tot.chi / rda2$tot.chi
  
  # Model SPEC
  spec_model <- rda(as.matrix(explanatory_merge[to_keep,climate_var]),explanatory_merge[to_keep,c("Long","Lat")])
  espc_model <- rda(explanatory_merge[to_keep,c("Long","Lat")],as.matrix(explanatory_merge[to_keep,climate_var]))
  
  
  # Return all elements
  out <- data.frame(climate_var=climate_var,
                    GSEC=part1$part$fract$R.squared[1],
                    GSEC_adj=part1$part$fract$Adj.R.squared[1],
                    GSSC=part1$part$fract$R.squared[2],
                    GSSC_adj=part1$part$fract$Adj.R.squared[2],
                    a=a,
                    b=ab - a,
                    c=c,
                    a_adj=part1$part$indfract$Adj.R.squared[1],
                    b_adj=part1$part$indfract$Adj.R.squared[2],
                    c_adj=part1$part$indfract$Adj.R.squared[3],
                    SPEC=RsquareAdj(spec_model)$r.squared,
                    SPEC_adj=RsquareAdj(spec_model)$adj.r.squared,
                    ESPC=RsquareAdj(espc_model)$r.squared,
                    ESPC_adj=RsquareAdj(espc_model)$adj.r.squared)
}

# If we removed any, add in a row of NAs so that we know...
if(length(clim_remove) > 0){
  for(i in 1:length(clim_remove)){
    variance_partitioning_clean <- rbind(variance_partitioning_clean,
                                         c(clim_remove[i],rep(NA,ncol(variance_partitioning_clean) - 1)))
  }
}


# Variance partitioning of all climate variables --------------------------
# This follows the workflow of Capblancq and Forester
message(">>> Starting all-climate variable modelling
        ")

## Null model
RDA0 <- rda(AF_mat_trans[to_keep,] ~ 1,  explanatory_merge[to_keep,]) 

## Full model
climate_vars = na.omit(variance_partitioning_clean)[,"climate_var"]
RDAfull <- rda(AF_mat_trans[to_keep,] ~ ., explanatory_merge[to_keep,c(climate_vars,"Lat","Long")])

## We need to check whether this model is actually able to be tested, and otherwise we'll build a new one...
if(is.na(RsquareAdj(RDAfull)$adj.r.squared)){
  # If true, the model has too many terms to be fitted.
  # To get around this, we'll fit terms sequentially based on largest > smallest GSEC to make the most complicated model we can
  climate_vars = variance_partitioning_clean[,"climate_var"][order(-variance_partitioning_clean[,"GSEC"])]
  model_vars = c(climate_vars[1],"Lat","Long")
  end_model = FALSE
  while(!end_model){
    RDAfull <- rda(AF_mat_trans[to_keep,],explanatory_merge[to_keep,model_vars])
    if(!is.na(RsquareAdj(RDAfull)$adj.r.squared)){
      # Add on one more climate var
      model_vars = c(model_vars[!model_vars %in% c("Lat","Long")],
                     climate_vars[!climate_vars %in% model_vars][1],
                     "Lat","Long")
    } else {
      end_model = TRUE
      final_model_vars = c(model_vars[1:(length(model_vars) - 3)],"Lat","Long")
      print(paste0("FINAL MODEL VARS: ",final_model_vars))
      
      # Fit the final model using formula notation
      reduced_explanatory_merge = explanatory_merge[,final_model_vars]
      RDAfull <- rda(AF_mat_trans[to_keep,] ~ ., data = reduced_explanatory_merge)
    }
  }
}

## Stepwise procedure with ordiR2step function
mod <- ordiR2step(RDA0, RDAfull, Pin = 0.05, R2permutations = 1000, R2scope = T)

# Save everything to an object
message(">>> Saving results and finishing...")
saveRDS(list(varpart = variance_partitioning_clean,
             full_model = RDAfull,
             reduced_model = mod,
             SNPN = nrow(AF_mat_impute)),
        paste0(res_dir,"/variance_partitioning_RDA_results.rds"))

