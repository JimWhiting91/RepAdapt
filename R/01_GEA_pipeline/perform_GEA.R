####################################################################################
# This script takes a VCF input and performs GEA analysis over climate variables
# Load these
lib <- c("hms","doParallel","parallel","data.table","sp","VGAM","argparse","raster","psych","pbmcapply","R.utils","tidyverse")
dummy <- suppressPackageStartupMessages(lapply(lib,function(x){ 
  if (!require(x,character.only = T)){ 
    install.packages(x,repos='https://utstat.toronto.edu/cran/')
    library(x,character.only = T)
  } else {
    library(x,character.only = T)
  }
}))

##########################################################################################
# Set up our arguments to parse
args <- commandArgs(asValues = T,excludeReserved = T)[-1]
print(args)

vcf_path <- args$vcf
n_cores <- args$n_cores
is.pool <- args$pool
output_dir <- args$dataset_dir

if(any(is.null(c(vcf_path,output_dir)))){
  stop("ERROR: No VCF/OUTPUT provided")
}
if(is.null(n_cores)){
  n_cores <- 1
}
if(is.null(is.pool)){
  is.pool <- FALSE
}

###########################################################################################
# # Test data... comment out as appropriate
# vcf_path <- "data/VCFs/19a_Ealbens_Murray/murray_Ealb_full_concatened_no_hybrids.vcf.gz"
# n_cores <- 12
# is.pool <- FALSE
# output_dir <- "220927_Eucalyptus_albens_Murray_Individual"

###########################################################################################

# Set up parallel cores
cl <- makeCluster(n_cores, type="FORK")
registerDoParallel(cl)
# Kill cluster
on.exit(stopCluster(cl))


###########################################################################################
# Report
message(paste0(">>> VCF = ",vcf_path,"
               "))
if(is.pool){
  message(">>> Treating data as PoolSeq
          ")
}

################################################################################################
# Function library

# Make a data.frame for climate data given lat long input and raster stack
make_climate_dataframe <- function(long_lat,climate_data){
  
  # Fetch points
  points <- SpatialPoints(long_lat, proj4string = climate_data@crs)
  values <- raster::extract(climate_data,points)
  
  # Pull data together
  df <- cbind.data.frame(coordinates(points),values)
  
  # Return the output
  return(df)
}

################################################################################################

# Based on the VCF, read in the metadata and subset it if needs be...
if(file.exists(paste0(dirname(vcf_path),"/sampling_data.csv"))){
  metadata <- read.csv(paste0(dirname(vcf_path),"/sampling_data.csv"))
} else {
  stop("Error: Expecting a file with sampling info as 'sampling_data.csv' in the same directory as the VCF, but this does not exist...")
}

# Make this if we haven't already...
dir.create("outputs/GEA_res",showWarnings = F)
# And make our outdir
dir.create(paste0("outputs/GEA_res/",output_dir),showWarnings = F)


# Fetch VCF individuals and further subset...
if(is.pool){
  # Pull column headers and clean
  vcf_inds <- system(paste0("head -n1 ",vcf_path),intern = T)
  vcf_inds <- unlist(strsplit(vcf_inds,"\t"))
  vcf_inds <- vcf_inds[grep(".FREQ",vcf_inds)]
  vcf_inds <- unique(gsub(".FREQ","",vcf_inds))
} else {
  vcf_inds <- system(paste0("bin/bcftools query -l ",vcf_path),intern = T)
}

# Do we need to clean?
if(any(!(vcf_inds %in% metadata$Sample_name))){
  to_remove <- names(table(c(vcf_inds,metadata$Sample_name)))[table(c(vcf_inds,metadata$Sample_name)) == 1]
  message(paste0(">>> Warning: The following individuals/pools are missing from VCF or sampling data: ",paste(to_remove,collapse=",")))
  
  metadata <- metadata[!(metadata$Sample_name %in% to_remove),]
  vcf_inds <- vcf_inds[!(vcf_inds %in% to_remove)]
  
}

# List final vcf inds
vcf_inds_final <- vcf_inds

# Report on this
species_code <- output_dir
workdir <- paste0("outputs/GEA_res/",species_code)

# Report start
message(">>> Starting analysis for ",species_code,"
")

# Remove any rows without lat long
metadata <- metadata[!(is.na(metadata$Lat)),]

# Group the individuals into populations based on geo-coords
metadata$geo_pops <- paste0(metadata$Lat,"_",metadata$Long)

# IF the data is pooled, check that we don't have duplicated geopops and if we do, bump them...
if(is.pool & any(duplicated(metadata$geo_pops))){
  message(">>> WARNING: Duplicated geopops detected with poolseq data. Only taking the first population for each geo-coord.
          ")
  pools_to_remove <- duplicated(metadata$geo_pops)
  metadata <- metadata[!(pools_to_remove),]
  
  # List final vcf inds
  vcf_inds_final <- vcf_inds[vcf_inds %in% metadata$Sample_name]
}

# Report
message(paste0(">>> ",species_code," contains ",length(unique(metadata$geo_pops))," geopops/pools
               "))

# For each of these geo-coord pops, subset the VCF and calculate allele frequencies
message(">>> Calculating Population AFs for Geo-Coords
        ")

# If data is based on individuals...
if(!file.exists(paste0(workdir,"/pop_allele_frqs.rds"))){
  if(!is.pool){
    geo_pop_freqs <- foreach(geo_pop = unique(metadata$geo_pops)) %dopar% {
      
      # Get the individuals we're getting frequencies for
      inds <- metadata[metadata$geo_pops == geo_pop,"Sample_name"]
      inds_to_keep <- paste0("--indv ",inds,collapse = " ")
      
      # Use VCFtools for frequency estimates
      system(paste0("vcftools --gzvcf ",vcf_path," --freq ",inds_to_keep," --out ",workdir,"/",geo_pop), ignore.stdout = TRUE, ignore.stderr = TRUE,wait = T)
      
      # Read back in the frequencies to R
      freq_res <- suppressMessages(data.frame(fread(paste0(workdir,"/",geo_pop,".frq"))))
      colnames(freq_res) <- c("chr","bp","allele_N","allele_count","ref_freq","alt_freq")
      
      # Remove this file
      system(paste0("rm -f ",workdir,"/",geo_pop,".frq"))
      
      # Reformat frequencies
      freq_res$ref_freq <- gsub("[^0-9.-]", "", freq_res$ref_freq)
      freq_res$alt_freq <- gsub("[^0-9.-]", "", freq_res$alt_freq)
      
      # Format for output
      freq_out <- data.frame(chr=freq_res$chr,
                             bp=freq_res$bp,
                             ref_freq=as.numeric(freq_res$ref_freq),
                             alt_freq=as.numeric(freq_res$alt_freq),
                             pop=geo_pop)
      return(freq_out)
    }
    
    # Save this list of frequencies to an RDS and remove the .frq files...
    saveRDS(geo_pop_freqs,paste0(workdir,"/pop_allele_frqs.rds"))
    
  } else if(is.pool){
    
    # Read in
    all_pool_AFs <- data.frame(fread(vcf_path,select=c("CHROM","POS",paste0(vcf_inds_final,".FREQ")),header=T))
    colnames(all_pool_AFs) <- c("chr","pos",vcf_inds_final)
    
    # Remove any bastard commas
    all_pool_AFs[,vcf_inds_final] <- apply(all_pool_AFs[,vcf_inds_final],2,gsub,pattern=",",replacement=".")
    
    # Sort these and organise to match the above
    geo_pop_freqs <- foreach(pool = metadata$Sample_name) %dopar% {
      
      # Pull from the whole
      freq_res <- all_pool_AFs[,c("chr","pos",pool)]
      colnames(freq_res)[3] <- "ref_freq"
      
      # Reformat frequencies
      freq_res$ref_freq <- as.numeric(gsub("%", "", freq_res$ref_freq))/100
      freq_res$alt_freq <- 1-freq_res$ref_freq
      
      # Make the geopop for this pool
      geopop <- paste0(metadata[metadata$Sample_name == pool,"Lat"],"_",metadata[metadata$Sample_name == pool,"Long"])
      
      # Format for output
      freq_out <- data.frame(chr=freq_res$chr,
                             bp=freq_res$pos,
                             ref_freq=as.numeric(freq_res$ref_freq),
                             alt_freq=as.numeric(freq_res$alt_freq),
                             pop=geopop)
      return(freq_out)
    }
    
    # Save these
    saveRDS(geo_pop_freqs,paste0(workdir,"/pop_allele_frqs.rds"))
    
  }
} else {
  geo_pop_freqs <- readRDS(paste0(workdir,"/pop_allele_frqs.rds"))
}

###### SUMMARISE POPULATIONS #######

###### BUILD CLINES #######
message(">>> Assembling Climate Clines Based on WorldClim 2.5")

# Fetch the climate data
dir.create("data/climate_data",showWarnings = F)
dir.create("data/climate_data/wordclim",showWarnings = F)
clima_stack <- raster::getData(name = 'worldclim', var = 'bio', res = 2.5, path = "data/climate_data/wordclim")

# Label the stack
names(clima_stack) <- c("mean_temp","mean_diurnal","isothermality","temp_seasonality","max_temp_warmest_month",
                        "min_temp_coldest_month","temp_range","mean_temp_wet_quarter","mean_temp_dry_quarter","mean_temp_warm_quarter",
                        "mean_temp_cold_quarter","annual_precip","precip_wet_month","precip_dry_month","precip_seasonality","precip_wet_quarter",
                        "precip_dry_quarter","precip_warm_quarter","precip_cold_quarter")

# Fetch our lat long co-ords again
# This maintains the order of long lat so they match the order of AF files imported
sub_long_lat <- data.frame(Long=as.numeric(sapply(strsplit(unique(metadata$geo_pops),"_"),'[[',2)),
                           Lat=as.numeric(sapply(strsplit(unique(metadata$geo_pops),"_"),'[[',1)))

# Make cline
species_cline <- make_climate_dataframe(sub_long_lat,clima_stack)

# Add to this estimates of climate change...
if(file.exists(paste0(dirname(vcf_path),"/climate_change_env.txt"))){
  climate_change = read.table(paste0(dirname(vcf_path),"/climate_change_env.txt"),header = T)
  climate_change = tidyr::separate(climate_change,"pop",sep = "_",into = c("Long","Lat"))
  species_cline = merge(species_cline,climate_change)
  
  # Don't forget to re-order
  species_cline = species_cline[order(match(paste0(species_cline$Lat,"_",species_cline$Long),
                                            paste0(sub_long_lat$Lat,"_",sub_long_lat$Long))),]
}

# Save the species cline
write.table(species_cline,paste0(workdir,"/climate_cline.tsv"),quote=F,row.names=F,sep="\t")

###### Associate AFs to climate #######
message(">>> Associating AFs to climate clines
        ")

# Build the AF matrix
af_mat <- matrix(ncol=nrow(geo_pop_freqs[[1]]),nrow=length(geo_pop_freqs))
for(i in 1:nrow(af_mat)){
  af_mat[i,] <- geo_pop_freqs[[i]][,"ref_freq"]
}
class(af_mat) <- "numeric"
rownames(af_mat) <- unique(metadata$geo_pops)
colnames(af_mat) <- paste0(geo_pop_freqs[[1]]$chr,":",geo_pop_freqs[[1]]$bp)

# Now make cor mats for all climate variables
climate_vars = colnames(species_cline)[!colnames(species_cline) %in% c("Long","Lat")]
for(x in 1:length(climate_vars)){
  
  # Message
  start_time <- Sys.time()
  message(paste0(">>> Starting GEA for ",climate_vars[x],"
                 "))
  
  # Set up
  tmp_var <- climate_vars[x]
  # Add a check here where if the GEA results already exist, and have the right number of SNPs, don't calculate them again...
  # Does the file exist?
  gea_exists = file.exists(paste0(workdir,"/",tmp_var,"_GEA.rds"))
  if(gea_exists){
    all_snps_check = nrow(readRDS(paste0(workdir,"/",tmp_var,"_GEA.rds"))) == ncol(af_mat)
  } else {
    all_snps_check = FALSE
  }
  
  # Now use these checks to work out what we have to run again...
  if(gea_exists & all_snps_check){
    message(">>> WARNING: Already GEA for all SNPs, skipping...")
  } else {
    if(gea_exists & !all_snps_check){
      message(">>> WARNING: File exists, but missing SNPs. Recalculating...")
    }
    
    cline_mat <- matrix(ncol=1,nrow=nrow(species_cline))
    cline_mat[,1] <- species_cline[,tmp_var]
    rownames(cline_mat) <- paste0(species_cline$Lat,"_",paste0(species_cline$Long))
    colnames(cline_mat) <- tmp_var
    
    # Reorder to match freqs
    if(any(rownames(cline_mat) != rownames(af_mat))){
      pops_to_keep = Reduce(intersect,list(rownames(af_mat),rownames(cline_mat)))
    } else {
      pops_to_keep = rownames(af_mat)
    }
    
    # Correlations in intervals
    interval=100
    winds=seq(1,ncol(af_mat),interval)
    winds2=winds+interval-1
    winds2[length(winds2)] <- ncol(af_mat)
    
    # Parallelise correlation calculations
    # cor_res <- foreach(y=1:length(winds)) %dopar% {
    cor_res <- pbmclapply(1:length(winds),function(y){
      
      # Psych method
      tmp_res <- psych::corr.test(x=cline_mat[pops_to_keep,1,drop=FALSE],
                                  y=af_mat[pops_to_keep,winds[y]:winds2[y]],
                                  method = "kendall",adjust="none")
      
      # Save an output for each climate var...
      tmp_out <- matrix(ncol=3,nrow=ncol(af_mat[,winds[y]:winds2[y]]))
      tmp_out[,1] <- tmp_res$r
      tmp_out[,2] <- tmp_res$p
      tmp_out[,3] <- tmp_res$se
      rownames(tmp_out) <- colnames(af_mat[,winds[y]:winds2[y]])
      colnames(tmp_out) <- c("tau.corr","p","se")
      return(tmp_out)
    },mc.cores=n_cores)
    # }
    
    # Bind these all together
    out <- data.frame(do.call(rbind,cor_res))
    out$snp_id <- rownames(out)
    out$climate_var <- tmp_var
    
    # Save these to the output dir
    saveRDS(out,paste0(workdir,"/",tmp_var,"_GEA.rds"))
    
    end_time <- Sys.time()
    message(paste0(">>> Finished GEA for ",climate_vars[x]," in ",as_hms(end_time-start_time)," H:M:S
                 "))
    
  }
  
}

# Send completion
message(">>> All Finished!
        ")


