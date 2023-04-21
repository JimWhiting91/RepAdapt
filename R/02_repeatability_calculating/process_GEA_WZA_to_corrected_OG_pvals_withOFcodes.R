# This script takes the results from GEA scans, and condenses them all down into a single table that shows for a given OG, Dataset, Climate Var, a single pval
lib <- c("poolr","regioneR","Rfast","mvmeta","qvalue","tidyr","ggtree","ape","VGAM","ggExtra","WeMix","pbmcapply","parallel","wCorr","data.table","ggplot2","viridis","ggridges","dplyr","patchwork","readr")
sapply(lib,library,character.only=T)
source("R/repadapt_functions.R")

############################################################################
# What GEA results are we looking at and processing...
run_name = "230321"
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"

# Where are our recombination maps
rec_maps = readRDS("outputs/processed_recombination_rates.rds")

############################################################################
# Assemble the Orthofinder Outputs to Analyse
n_cores = 6 # Number of cpu
orthogroup_cutoff <- 10 # Maximum number of paralogs in an OG
snp_bin_size <- 100 # How many genes in each bin for estimating wza sd
min_snp_per_gene <- 5

# Where are the blast results?
OG_dir <- "outputs/orthology/Results_221213_18_genomes_Ptaeda_isoforms_removed/Orthogroups/"

# Where is all the metadata...
vcf_genome_map <- read.table("metadata/vcf_genome_gff_220830_map.txt",fill=T)
dataset_meta <- vcf_genome_map[,c("V8","V5")]
colnames(dataset_meta) <- c("dataset","genome")

genome_OF_codes <- data.frame(prot=c("Aalpina.faa",
                                     "Ahalleri.faa",
                                     "Alyrata.faa",
                                     "Athaliana.faa",
                                     "Atubercatus.faa",
                                     "Bstricta.faa",
                                     "Crubella.faa",
                                     "Egrandis.faa",
                                     "Hannuus.faa",
                                     "Mtruncatula.faa",
                                     "Pabies.faa",
                                     "Pdeltoides.faa",
                                     "Phallii.faa",
                                     "Pmenziesii.faa",
                                     "Ptaeda.faa",
                                     "Ptremula.faa",
                                     "Ptrichocarpa.faa",
                                     "Qpetraea.faa"))
genome_OF_codes$genome <- gsub(".faa","",genome_OF_codes$prot)
genome_OF_codes$OF_code <- (1:nrow(genome_OF_codes))-1

# Read in all of our liftovers from gff to OF codes
gff_OF_liftovers <- list.files("data/reference_genomes/",recursive = TRUE,pattern="_proteome_to_OF_id_map.txt")
gff_OF_liftovers <- lapply(genome_OF_codes$genome,function(x) grep(x,gff_OF_liftovers,value=T))
names(gff_OF_liftovers) <- genome_OF_codes$genome

############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
############################################################################################################################################################
# Prepare Orthofinder data

# Fetch the sequence IDs...
OG_sequences <- read.table(paste0(dirname(OG_dir),"/WorkingDirectory/SequenceIDs.txt"),fill=T)[,1:2]

# Remove empties
OG_sequences <- OG_sequences[OG_sequences$V2 != "",]

# Format text
colnames(OG_sequences) <- c("OF_gene","gene_name")
OG_sequences$OF_gene <- gsub(":","",OG_sequences$OF_gene)
OG_sequences$OF_genome <- sapply(strsplit(OG_sequences$OF_gene,"_"),'[[',1)
OG_sequences$gene_name_noGenome <- OG_sequences$gene_name
for(genome in genome_OF_codes$genome){
  OF_code <- genome_OF_codes[genome_OF_codes$genome == genome,"OF_code"]
  OG_sequences[OG_sequences$OF_genome == as.character(OF_code),"gene_name"] <- paste0(genome,"_",OG_sequences[OG_sequences$OF_genome == as.character(OF_code),"gene_name"])
}

# Transform names so they match
OG_sequences$gene_name <- gsub("transcript:","transcript_",OG_sequences$gene_name)
OG_sequences$gene_name_noGenome <- gsub("transcript:","transcript_",OG_sequences$gene_name_noGenome)

OG_sequences$gene_name <- gsub("RNA:","RNA_",OG_sequences$gene_name)
OG_sequences$gene_name_noGenome <- gsub("RNA:","RNA_",OG_sequences$gene_name_noGenome)

# We also now want to attach the Orthogroup to each OG_sequence...
all_OG <- data.frame(fread(paste0(OG_dir,"Orthogroups.tsv")))
all_OG_long <- rbindlist(pbmclapply(2:ncol(all_OG),function(col){

  # Subset...
  tmp <- data.frame(separate_rows(all_OG[,c(1,col)],colnames(all_OG)[col],sep=","))
  tmp[,2] <- gsub(" ","",tmp[,2])
  tmp$genome <- colnames(tmp)[2]
  colnames(tmp)[2] <- "gene_name_noGenome"
  tmp$gene_name <- paste0(tmp$genome,"_",tmp$gene_name_noGenome)
  return(tmp[tmp$gene_name != paste0(tmp$genome[1],"_"),])
},mc.cores = n_cores))

# Merge it with OG_sequences
OG_sequences2 <- merge(OG_sequences,all_OG_long[,c("Orthogroup","gene_name")],by="gene_name")
OG_sequences2$full_OF <- OG_sequences2$OF_gene

for(OF_code in genome_OF_codes$OF_code){
  OG_sequences2[OG_sequences2$OF_genome == OF_code,"genome"] <- genome_OF_codes[genome_OF_codes$OF_code == OF_code,"genome"]
}

# Add in the original gea co-ordinates
all_gff_map = rbindlist(lapply(names(gff_OF_liftovers),function(genome_tmp){

  # Read in gff map
  gff_tmp = fread(paste0("data/reference_genomes/",gff_OF_liftovers[genome_tmp]))

  # Edit the OF_ID to add the genome code
  gff_tmp$full_OF = paste0(genome_OF_codes$OF_code[which(genome_OF_codes$genome == genome_tmp)],"_",gff_tmp$OF_ID)

  colnames(gff_tmp)[3] = "seqname"
  return(gff_tmp)
}))

# Merge with the OG_sequences...
gea_gene_OF_map = merge(all_gff_map,OG_sequences2[,c("full_OF","genome","Orthogroup")],by = "full_OF",all.x = T)
gea_gene_OF_map$gea_gene = paste0(gea_gene_OF_map$seqname,":",gea_gene_OF_map$start,"-",gea_gene_OF_map$end)

# Save this map
saveRDS(gea_gene_OF_map,
            paste0("data/OF_OG_gea_gene_map_",output_name,".rds"))

############################################################################################################################################################
# Prepare GEA data -------------------------------------------------------
focal_datasets <- list.files("outputs/GEA_res",pattern = run_name)
focal_datasets <- grep("processed",focal_datasets,invert=T,value=T)
focal_datasets <- grep("CoAdapTree",focal_datasets,invert = T,value=T)
focal_datasets <- grep(".rds",focal_datasets,invert = T,value=T)

############################################################################################################################################################
# Define climate variaibles here
climate_vars <- c("mean_temp","mean_diurnal","isothermality","temp_seasonality","max_temp_warmest_month",
                  "min_temp_coldest_month","temp_range","mean_temp_wet_quarter","mean_temp_dry_quarter","mean_temp_warm_quarter",
                  "mean_temp_cold_quarter","annual_precip","precip_wet_month","precip_dry_month","precip_seasonality","precip_wet_quarter",
                  "precip_dry_quarter","precip_warm_quarter","precip_cold_quarter","tmax_clim_change","prec_clim_change")

all_climate_res <- lapply(climate_vars,function(focal_climate){
  
  message(paste0("STARTING ",focal_climate))

  # First remove any datasets where we don't have GEA results...
  focal_datasets_climate <- na.omit(sapply(focal_datasets,function(x){
    if(file.exists(paste0("outputs/GEA_res/",x,"/",focal_climate,"_WZA_TC_allgenes.rds"))){
      return(x)
    } else {
      return(NA)
    }
  }))

  # Read in all of our GEA results
  message(">>> Fetching GEA WZA scores")

  # Loop over datasets
  all_gea_res <- lapply(focal_datasets_climate,function(dataset){

    print(dataset)

    # Get raw GEA
    wza <- readRDS(paste0("outputs/GEA_res/",dataset,"/",focal_climate,"_WZA_TC_allgenes.rds"))

    # Merge these with the relevant OF_code...
    genome = na.omit(dataset_meta[dataset_meta$dataset == gsub(paste0(run_name,"_"),"",dataset),"genome"])
    wza$full_OF = paste0(which(genome == unique(sort(dataset_meta$genome))) - 1,"_",wza$OF_ID)
    
    # Fetch the orthogroups as well...
    wza_merge <- merge(wza,OG_sequences2[,c("Orthogroup","full_OF")],"full_OF")
    OG_counts1 <- table(wza_merge$Orthogroup)
    wza_merge <- wza_merge[Orthogroup %in% names(OG_counts1[OG_counts1 <= orthogroup_cutoff]),]

    # Calculate parametric pvalues by sd modelling...
    snp_groups <- rep(1:floor(nrow(wza_merge)/snp_bin_size),each=snp_bin_size)
    to_model <- wza_merge[order(snp_count),][1:length(snp_groups),]
    to_model$snp_group <- snp_groups

    #### Do the mirroring ####
    snp_bin_sd <- to_model[,mirror_wza_sd_calc(wza_vector = .SD$weiZ_downsample,
                                               snp_count_vector = .SD$downsampled_snp_count),
                           by=snp_group]

    #### Continue ####
    snp_bin_mean <- to_model[,.(snp_count = mean(downsampled_snp_count),
                                mean_wza = mean(weiZ_downsample)),by=snp_group]

    # Fit a spline to this....
    mean_spline = smooth.spline(x = log10(snp_bin_mean$snp_count),
                                y = snp_bin_mean$mean_wza,spar = 1)
    sd_spline = smooth.spline(x = log10(snp_bin_sd$snp_count),
                              y = snp_bin_sd$wza_sd,spar = 1)

    # Use this to predict sd and pnorm some pvals
    wza_merge$sd_pvalue <- pnorm(wza_merge$weiZ_downsample,
                                 # mean = mean(wza_merge$weiZ_downsample,na.rm=T),
                                 mean = predict(log10(wza_merge$downsampled_snp_count),object=mean_spline)$y,
                                 sd = predict(log10(wza_merge$downsampled_snp_count),object=sd_spline)$y,
                                 lower.tail = F)

    # Mark
    wza_merge$dataset <- gsub(paste0(run_name,"_"),"",dataset)
    wza_merge$genome <- genome

    # Now remove the smallest genes that are < cutoff
    if(min_snp_per_gene >= quantile(wza_merge$snp_count,0.05)){
      wza_merge <- wza_merge[wza_merge$snp_count >= quantile(wza_merge$snp_count,0.05),]
    } else {
      wza_merge <- wza_merge[wza_merge$snp_count >= min_snp_per_gene,]
    }

    # Add an empirical p-value for SNP count...
    wza_merge$snp_count_pval <- rank(-wza_merge$downsampled_snp_count,na.last = "keep")/(sum(is.na(wza_merge$downsampled_snp_count) == F)+1)

    # Finally convert the sd pvalue to the sd epvalue
    wza_merge$sd_epvalue <- rank(wza_merge$sd_pvalue)/(sum(!is.na(wza_merge$sd_pvalue))+1)

    # And get the species name...
    wza_merge$species = paste(strsplit(wza_merge$dataset[1],"_")[[1]][1:2],collapse = " ")


    ############################ Add recombination rates if we can... ############################
    if(genome %in% names(rec_maps)){
      rec_tmp = rec_maps[[genome]]

      # Intersect...
      if(any(colnames(wza_merge) == 'gene_id')){
        wza_merge$gea_gene = wza_merge$gene_id
      } 
      # Gene regions
      gene_regions = tidyr::separate(data.frame(wza_merge[,"gea_gene"]), col = "gea_gene",
                                     into = c("chr", "start", "end"), sep = ":|-")
      gene_regions$start = as.integer(gene_regions$start) - 500
      gene_regions = toGRanges(gene_regions)

      # Rec regions
      rec_tmp$rec_region = paste0(rec_tmp$chr,":",rec_tmp$start,"-",rec_tmp$end)
      rec_regions = toGRanges(rec_tmp[,c("chr","start","end")])

      # Overlap
      rec_gene_overlap = overlapRegions(gene_regions,rec_regions,get.bases = T)

      # Merge with rec_info...
      rec_gene_overlap$gea_gene = paste0(rec_gene_overlap$chr,":",rec_gene_overlap$startA + 500,"-",rec_gene_overlap$endA)
      rec_gene_overlap$rec_region = paste0(rec_gene_overlap$chr,":",rec_gene_overlap$startB,"-",rec_gene_overlap$endB)
      rec_gene_overlap_merge = merge(rec_gene_overlap,rec_tmp[,.(rec_region,cM_Mb)],by="rec_region")

      # Average through within gea_gene
      rec_gene_rates = data.table(rec_gene_overlap_merge)[,.(rec_rate = weighted.mean(.SD$cM_Mb,w = .SD$ov.bases/sum(.SD$ov.bases))),by=.(gea_gene)]

      # Add these back onto wza_merge
      wza_merge = data.frame(merge(wza_merge,rec_gene_rates,by="gea_gene",all.x=T))
      wza_merge$rec_rate_pval = NA

      # Set empirical pvalues explicitly such that low rec rates = low pvals
      wza_merge[!is.na(wza_merge$rec_rate),"rec_rate_pval"] = empPvals( -wza_merge[!is.na(wza_merge$rec_rate),"rec_rate"],
                                                                        -wza_merge[!is.na(wza_merge$rec_rate),"rec_rate"])

      ##############################################################################################

    } else {
      wza_merge$rec_rate = NA
      wza_merge$rec_rate_pval = NA
    }

    # Estimate rec rate empirical pvals...


    return(data.table(wza_merge[,c("species","dataset","Orthogroup","mean_corr","sd_pvalue","sd_epvalue","weiZ","weiZ_downsample","full_OF","snp_count","downsampled_snp_count","snp_count_pval","rec_rate_pval")]))
  })
  names(all_gea_res) <- focal_datasets_climate
  
  # Save these at this stage so we can do the recombination rate bits...
  saveRDS(na.omit(rbindlist(all_gea_res)),
          paste0("outputs/GEA_res/run",run_name,"_",output_name,"_RecRate_",focal_climate,"_pvals.rds"))


  #### Here we combine pvals for the same species for the same gene ####
  # BUT DON'T DO THE EMPIRICAL STEP
  OG_combinedpvals = unique(rbindlist(all_gea_res)[,.(combined_sd_pvalue = poolr::fisher(.SD$sd_pvalue)$p,
                                                      Orthogroup), by = .(full_OF,species)])
  OG_combinedpvals$climate_var = focal_climate

  return(OG_combinedpvals)
})
names(all_climate_res) = climate_vars

# Within Orthogroups and Datasets, shuffle climates...
saveRDS(all_climate_res,
        paste0("outputs/GEA_res/run",run_name,"_",output_name,"_WZA_OG_PerGene_pvals.rds"))
