# Fetch specific gene info driving repeatability
lib <- c("biomaRt","org.At.tair.db","ggridges","ggrepel","cowplot","ggplot2","data.table","tidyr","pbmcapply","ggtree","ape","viridis","STRINGdb","qvalue","dplyr")
sapply(lib,library,character.only=T)
n_cores = 6

# Fetch drivers
run_name = "230321"
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"
drivers = readRDS(paste0('outputs/',output_name,'_dataset_climate_signif.rds'))

# And the original all gene p-values
pvals_file = readRDS(paste0("outputs/GEA_res/run",run_name,"_",output_name,"_WZA_OG_PerGene_pvals.rds"))

# We'll also need all of the OF maps...
OF_map = readRDS(paste0("data/OF_OG_gea_gene_map_",output_name,".rds"))

# Now just loop through the drivers and return the relevant gene...
driving_genes = lapply(1:nrow(drivers),function(x){
  print(x)
  
  OG_tmp = drivers$Orthogroup[x]
  clim_tmp = drivers$climate_var[x]
  species_tmp = drivers$dataset[x]
  
  # Fetch the specific gene and its original pval
  tmp = pvals_file[[clim_tmp]][species == species_tmp & Orthogroup == OG_tmp,][order(combined_sd_pvalue),][1,]
  # Attach the genome info
  cbind(tmp,
        OF_map[full_OF == tmp$full_OF[1],.(seqname,start,end,gea_gene,genome)])
}) |>
  rbindlist()

# Save this
write.csv(driving_genes,
          paste0('outputs/',output_name,'_specific_genes_driving_repeatability.csv'),
          quote = F,row.names = F)
