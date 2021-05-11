####################################################################################
# This script takes a VCF input and performs GEA analysis over climate variables
lib <- c("parallel","data.table","sp","VGAM")
sapply(lib,library,character.only=T)

# Read the VCF in from the command line
args <- commandArgs(TRUE)
vcf_path <- as.character(args[1])
metadata_path <- as.character(args[2])
n_cores <- as.character(args[3])

# For setting up purposes
vcf_path <- "/lu213/james.whiting/RepAdapt/data/VCFs/06_Ahalleri_Kubota_wgs/Ahalg_full_concatened.vcf.gz"
metadata_path <- "metadata/sample_species_vcf_author_map_v2_210507.csv"
n_cores <- 16
################################################
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

################################################

# Based on the VCF, read in the metadata and subset it
metadata <- read.csv(metadata_path)
sub_meta <- metadata[metadata$VCF==vcf_path,]

# Make a working dir
sub_meta$Taxon_name <- gsub(" ","_",sub_meta$Taxon_name)
dir.create(paste0("outputs/",sub_meta$Taxon_name[1],"_",sub_meta$Author[1],"_",sub_meta$Data_type[1]))
workdir <- paste0("outputs/",sub_meta$Taxon_name[1],"_",sub_meta$Author[1],"_",sub_meta$Data_type[1])

# Group the individuals into populations based on geo-coords
sub_meta$geo_pops <- paste0(sub_meta$Lat,"_",sub_meta$Long)

# For each of these geo-coord pops, subset the VCF and calculate allele frequencies
cat("Calculating Population AFs for Geo-Coords")
geo_pop_freqs <- mclapply(unique(sub_meta$geo_pops),function(geo_pop){

  # Get the individuals we're getting frequencies for
  inds <- sub_meta[sub_meta$geo_pops == geo_pop,"Sample_name"]
  inds_to_keep <- paste0("--indv ",inds,collapse = " ")

  # Use VCFtools for frequency estimates
  if(!(file.exists(paste0(workdir,"/",geo_pop,".frq")))){
    system(paste0("vcftools --gzvcf ",vcf_path," --freq ",inds_to_keep," --out ",workdir,"/",geo_pop),wait = T)
  }

  # Read back in the frequencies to R
  freq_res <- data.frame(fread(paste0(workdir,"/",geo_pop,".frq")))
  colnames(freq_res) <- c("chr","bp","allele_N","allele_count","ref_freq","alt_freq")

  # Format for output
  freq_out <- data.frame(chr=freq_res$chr,
                         bp=freq_res$bp,
                         ref_freq=sapply(strsplit(freq_res$ref_freq,":"),'[[',2),
                         alt_freq=sapply(strsplit(freq_res$alt_freq,":"),'[[',2),
                         pop=geo_pop)
  return(freq_out)
},mc.cores=n_cores)

###### SUMMARISE POPULATIONS #######

###### BUILD CLINES #######
cat("Assembling Climate Clines Based on WorldClim 2.5")

# Fetch the climate data
dir.create("data/climate_data")
dir.create("data/climate_data/wordclim")
clima_stack <- raster::getData(name = 'worldclim', var = 'bio', res = 2.5, path = "data/climate_data/wordclim")

# Label the stack
names(clima_stack) <- c("mean_temp","mean_diurnal","isothermality","temp_seasonality","max_temp_warmest_month",
                        "min_temp_coldest_month","temp_range","mean_temp_wet_quarter","mean_temp_dry_quarter","mean_temp_warm_quarter",
                        "mean_temp_cold_quarter","annual_precip","precip_wet_month","precip_dry_month","precip_seasonality","precip_wet_quarter",
                        "precip_dry_quarter","precip_warm_quarter","precip_cold_quarter")
climate_vars <- c("Annual Mean Temperature",
                  "Mean Diurnal Range",
                  "Isothermality",
                  "Temperature Seasonality",
                  "Max Temperature of Warmest Month",
                  "Min Temperature of Coldest Month",
                  "Temperature Annual Range",
                  "Mean Temperature of Wettest Quarter",
                  "Mean Temperature of Driest Quarter",
                  "Mean Temperature of Warmest Quarter",
                  "Mean Temperature of Coldest Quarter",
                  "Annual Precipitation",
                  "Precipitation of Wettest Month",
                  "Precipitation of Driest Month",
                  "Precipitation Seasonality",
                  "Precipitation of Wettest Quarter",
                  "Precipitation of Driest Quarter",
                  "Precipitation of Warmest Quarter",
                  "Precipitation of Coldest Quarter")

# Fetch our lat long co-ords again
# This maintains the order of long lat so they match the order of AF files imported
sub_long_lat <- data.frame(Long=as.numeric(sapply(strsplit(unique(sub_meta$geo_pops),"_"),'[[',2)),
                           Lat=as.numeric(sapply(strsplit(unique(sub_meta$geo_pops),"_"),'[[',1)))

# Make cline
species_cline <- make_climate_dataframe(sub_long_lat,clima_stack)

# Save the species cline
write.table(species_cline,paste0(workdir,"/climate_cline.tsv"),quote=F,row.names=F,sep="\t")

###### Associate AFs to climate #######
cat("Associating AFs to climate clines")

# Build the AF matrix
af_mat <- matrix(ncol=nrow(geo_pop_freqs[[1]]),nrow=length(geo_pop_freqs))
for(i in 1:nrow(af_mat)){
  af_mat[i,] <- geo_pop_freqs[[i]][,"ref_freq"]
}
class(af_mat) <- "numeric"
rownames(af_mat) <- unique(sub_meta$geo_pops)
colnames(af_mat) <- paste0(geo_pop_freqs[[1]]$chr,":",geo_pop_freqs[[1]]$bp)

# Now make cor mats for all climate variables
cline_corrs <- lapply(1:length(climate_vars),function(x){

  # Message
  message(paste0("Starting GEA for ",climate_vars[x]))

  # Set up
  tmp_var <- names(clima_stack)[x]
  cline_mat <- matrix(ncol=1,nrow=nrow(species_cline))
  cline_mat[,1] <- species_cline[,tmp_var]
  rownames(cline_mat) <- rownames(af_mat)
  colnames(cline_mat) <- tmp_var

  # Correlations in intervals
  interval=1000
  winds=seq(1,ncol(af_mat),interval)
  winds2=winds+interval-1
  winds2[length(winds2)] <- ncol(af_mat)
  cor_res <- pbmcapply::pbmclapply(1:length(winds),function(y){
    tmp_res <- psych::corr.test(x=cline_mat,y=af_mat[,winds[y]:winds2[y]],method = "kendall",adjust="none")
    tmp_out <- matrix(ncol=3,nrow=ncol(af_mat[,winds[y]:winds2[y]]))
    tmp_out[,1] <- tmp_res$r
    tmp_out[,2] <- tmp_res$p
    tmp_out[,3] <- tmp_res$se
    rownames(tmp_out) <- colnames(af_mat[,winds[y]:winds2[y]])
    colnames(tmp_out) <- c("tau.corr","p","se")
    return(tmp_out)
  },mc.cores=n_cores)

  # Bind these all together
  out <- data.frame(do.call(rbind,cor_res))
  out$snp <- rownames(out)
  out$climate_var <- tmp_var

  # Save these to the output dir
  write.table(out,
              paste0(workdir,"/",tmp_var,"_GEA.tsv"),
              sep="\t",quote=F)
})
