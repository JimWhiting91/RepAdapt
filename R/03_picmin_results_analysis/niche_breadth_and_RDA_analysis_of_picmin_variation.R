#### This script brings together RDA within each dataset to gauge which variables are most confounded by spatial autocorrelation
lib = c("scatterplot3d","cowplot","ggridges","viridis","tidyr","ape","DHARMa","sjPlot","qvalue","ggplot2","pbmcapply","biomaRt","org.At.tair.db","S4Vectors", "IRanges", "GenomicRanges", "SummarizedExperiment","readr","data.table")
sapply(lib,library,character.only=T)
n_cores = 6

# Which run name are we looking at?
run_name = "230321"
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"
picmin_outputs = readRDS(paste0("outputs/",output_name,"_picmin_results_doubleDS_pvals.rds"))
picmin_fdr = rbindlist(lapply(picmin_outputs,'[[',1))
OG_pvals = rbindlist(lapply(picmin_outputs,'[[',2))

# Extract and clean FDR results -------------------------------------------
climate_vars = unique(picmin_fdr$climate_var)

# Identify strongest repeatability signatures...
repeatable_OG_climate = picmin_fdr[picmin_p_adj < 0.005 & 
                                     !climate_var %in% c('tmax_clim_change','prec_clim_change'),]

# Run through species and count up the number of times they contribute to convergence signals. This is akin to their Fig 2B values.
# We want to count the per species-climate contributions
species_climate_contributions = rbindlist(pbmclapply(1:nrow(repeatable_OG_climate),function(x){
  pvals_tmp = OG_pvals[Orthogroup == repeatable_OG_climate$Orthogroup[x] &
                         climate_var == repeatable_OG_climate$climate[x],.(species,epval_final,climate_var,Orthogroup)][order(epval_final),]
  pval_cutoff = pvals_tmp$epval_final[repeatable_OG_climate$config_est[x]]
  if(pval_cutoff > 0.1){
    pvals_tmp = pvals_tmp[epval_final < 0.1,]
  } else {
    pvals_tmp = pvals_tmp[epval_final <= pval_cutoff,]
  }
  return(pvals_tmp)
},mc.cores = n_cores))

# How many Orthogroups is this per climate variable...
repeatable_OG_counts = unique(species_climate_contributions[,.(Orthogroup,climate_var)])[,.(OG_N = nrow(.SD)),by = climate_var][order(OG_N),]
# And how many is it per species?
repeatable_OG_counts_species = unique(species_climate_contributions[,.(Orthogroup,species)])[,.(species_OG_N = nrow(.SD)),by = species][order(species_OG_N),]

# Count up the contributions by species by climate var
species_climate_contributions = species_climate_contributions[,.(contributing_N = nrow(.SD)),by = .(climate_var,species)]
ggplot(species_climate_contributions,aes(contributing_N)) +
  geom_histogram() +
  facet_wrap(~species)

# Merge with counts and also estimate the proportion...
species_climate_contributions = merge(species_climate_contributions,repeatable_OG_counts)
species_climate_contributions = merge(species_climate_contributions,repeatable_OG_counts_species,by = "species")

species_climate_contributions$contributing_prop = species_climate_contributions$contributing_N/species_climate_contributions$OG_N
species_climate_contributions$contributing_prop_species = species_climate_contributions$contributing_N/species_climate_contributions$species_OG_N

# Fetch and prepare RDA outputs -------------------------------------------
# Find all of the RDA outputs...
res_dir = grep(".rds",list.files("outputs/GEA_res/",run_name),invert = T,value = T)
rda_res = rbindlist(lapply(paste0("outputs/GEA_res/",res_dir,"/variance_partitioning_RDA_results.rds"),function(x){
  tmp = data.frame(readRDS(x)$varpart)
  for(i in 2:ncol(tmp)){
    tmp[,i] = as.numeric(tmp[,i])
  }
  tmp$dataset = basename(dirname(x))
  tmp
}))
rda_res[SPEC_adj < 0, SPEC_adj := 0]

# Visualise the space, env, confounded proportions...
to_plot = melt(rda_res[,.(dataset,climate_var,a,b,c)])
dataset_GSSC = unique(na.omit(rda_res[,.(dataset,GSSC)]))
to_plot$dataset = gsub(paste0(run_name,"_"),"",to_plot$dataset)
dataset_GSSC$dataset = gsub(paste0(run_name,"_"),"",dataset_GSSC$dataset)
ggplot(to_plot,aes(y = climate_var,x = value,fill = variable)) +
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~dataset) +
  geom_vline(data = dataset_GSSC,aes(xintercept = GSSC))

# Build a species-dataset map
species_dataset_map = data.table(dataset = unique(rda_res$dataset))
species_dataset_map$species = gsub(paste0(run_name,"_"),"",species_dataset_map$dataset)
species_dataset_map$species = paste0(sapply(strsplit(species_dataset_map$species,"_"),'[[',1)," ",sapply(strsplit(species_dataset_map$species,"_"),'[[',2))

# For this analysis, we are primarily interested in GSEC, i.e. how much of the genome is associated with environment...
gsec_res = rda_res[,.(climate_var,dataset,GSEC)]
# Add species
gsec_res = merge(gsec_res,species_dataset_map)

# Compare RDA and picmin contributions... ---------------------------------
gsec_picmin_merge = merge(species_climate_contributions,gsec_res,by = c("species","climate_var"))

# In each species, what is the association between GSEC and contribution to picmin signals?
# Plot these
# gsec_picmin_merge$climate_var_F = stringr::str_to_title(gsub("_"," ",gsec_picmin_merge$climate_var))
GSEC_fig = ggplot(gsec_picmin_merge,aes(x = GSEC,y = contributing_prop)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~species,scales = "free",strip.position = "top") +
  theme_bw() +
  theme() +
  labs(x = "Proportion of population structure var explained by climate",
       y = "Proportional contributions to repeatability") 

# Prepare niche breadth for similar comparison ----------------------------
# Fetch the original climate_clines
climate_clines = rbindlist(lapply(paste0("outputs/GEA_res/",res_dir,"/climate_cline.tsv"),function(x){
  tmp = fread(x)
  tmp$dataset = basename(dirname(x))
  tmp
}))
climate_clines_long = melt(climate_clines[,c("dataset",climate_vars),with = FALSE])
colnames(climate_clines_long) = c("dataset","climate_var","value")

# Calculate global max and mins, with tail trimming as well...
global_max_mins = climate_clines_long[,.(global_min = min(value,na.rm = T),
                                         global_max = max(value,na.rm = T),
                                         global_minQ = quantile(na.omit(value),probs = 0.1),
                                         global_maxQ = quantile(na.omit(value),probs = 0.9)),by = climate_var]
global_max_mins$global_range = global_max_mins$global_max - global_max_mins$global_min
global_max_mins$global_rangeQ = global_max_mins$global_maxQ - global_max_mins$global_minQ

# Do the same for each individual variable...
dataset_clim_max_mins = climate_clines_long[,.(clim_min = min(value,na.rm = T),
                                               clim_max = max(value,na.rm = T),
                                               clim_minQ = quantile(na.omit(value),probs = 0.1),
                                               clim_maxQ = quantile(na.omit(value),probs = 0.9)),by = .(dataset,climate_var)]
dataset_clim_max_mins$diffs = dataset_clim_max_mins$clim_max - dataset_clim_max_mins$clim_min
dataset_clim_max_mins$diffsQ = dataset_clim_max_mins$clim_maxQ - dataset_clim_max_mins$clim_minQ

# Standardise on the basis of global values
dataset_clim_max_mins = merge(dataset_clim_max_mins,global_max_mins[,.(climate_var,global_range,global_rangeQ)],by = "climate_var")
dataset_clim_max_mins$relative_breadth = dataset_clim_max_mins$diffs/dataset_clim_max_mins$global_range
dataset_clim_max_mins$relative_breadthQ = dataset_clim_max_mins$diffsQ/dataset_clim_max_mins$global_rangeQ

# Similarly to above, merge these with the estimates of picmin contributions...
gsec_picmin_NB_merge = merge(gsec_picmin_merge,dataset_clim_max_mins[,.(dataset,climate_var,relative_breadth,relative_breadthQ)],by = c("dataset","climate_var"))

# Plot these
gsec_picmin_NB_merge$climate_var_F = stringr::str_to_title(gsub("_"," ",gsec_picmin_NB_merge$climate_var))
NB1_fig = ggplot(gsec_picmin_NB_merge,aes(x = relative_breadth,y = contributing_prop)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~climate_var_F,scales = "free",strip.position = "top") +
  theme_bw() +
  theme() +
  labs(x = "Niche Breadth (Relative to Global)",
       y = "Proportional contributions to repeatability") +
  ggtitle("Niche Breadth - Standardised range relative to global range")



# Repeat the above, but standardise by means instead ----------------------
# For proportional variables: we need to replace the mean with (p * 1-p)
# isothermality and precip_seasonality
proportional_vars = c("isothermality","precip_seasonality")

# Also transform temp vars to Kelvin...
temps_to_transform = c(grep("seasonality",grep("temp",climate_vars,value = T),invert = T,value = T),"mean_diurnal")
climate_clines_long_trans = climate_clines_long
climate_clines_long_trans[climate_var %in% temps_to_transform,value := value/10 +  273.15] 
ggplot(climate_clines_long_trans,aes(value)) + geom_histogram() + facet_wrap(~climate_var,scales = "free")
ggplot(climate_clines_long_trans[climate_var %in% temps_to_transform,],aes(value)) + geom_histogram() + facet_wrap(~climate_var,scales = "free")


# Fetch values
dataset_clim_max_mins = climate_clines_long_trans[!climate_var %in% proportional_vars,][,.(clim_min = min(value,na.rm = T),
                                                                                           clim_max = max(value,na.rm = T),
                                                                                           clim_minQ = quantile(na.omit(value),probs = 0.1),
                                                                                           clim_maxQ = quantile(na.omit(value),probs = 0.9),
                                                                                           clim_mean = mean(value,na.rm = T)),by = .(dataset,climate_var)]
dataset_clim_max_mins$diffs = (dataset_clim_max_mins$clim_max - dataset_clim_max_mins$clim_min) / dataset_clim_max_mins$clim_mean
dataset_clim_max_mins$diffsQ = (dataset_clim_max_mins$clim_maxQ - dataset_clim_max_mins$clim_minQ) / dataset_clim_max_mins$clim_mean

dataset_clim_max_mins_props = climate_clines_long[climate_var %in% proportional_vars,][,.(clim_min = 0.01 * min(value,na.rm = T),
                                                                                          clim_max = 0.01 * max(value,na.rm = T),
                                                                                          clim_minQ = 0.01 * quantile(na.omit(value),probs = 0.1),
                                                                                          clim_maxQ = 0.01 * quantile(na.omit(value),probs = 0.9)),by = .(dataset,climate_var)]

dataset_clim_max_mins_props$diffs = (dataset_clim_max_mins_props$clim_max - dataset_clim_max_mins_props$clim_min)
dataset_clim_max_mins_props$diffsQ = (dataset_clim_max_mins_props$clim_maxQ - dataset_clim_max_mins_props$clim_minQ)
dataset_clim_max_mins_props$diffs = dataset_clim_max_mins_props$diffs * (1 - dataset_clim_max_mins_props$diffs)
dataset_clim_max_mins_props$diffsQ = dataset_clim_max_mins_props$diffsQ * (1 - dataset_clim_max_mins_props$diffsQ)

# Bind back together...
dataset_clim_max_mins = rbind(dataset_clim_max_mins[,.(dataset,climate_var,diffs,diffsQ)],
                              dataset_clim_max_mins_props[,.(dataset,climate_var,diffs,diffsQ)])


# Similarly to above, merge these with the estimates of picmin contributions...
gsec_picmin_NB_merge = merge(gsec_picmin_merge,dataset_clim_max_mins[,.(dataset,climate_var,diffs,diffsQ)],by = c("dataset","climate_var"))

# Plot these
gsec_picmin_NB_merge$climate_var_F = stringr::str_to_title(gsub("_"," ",gsec_picmin_NB_merge$climate_var))
NB2_fig = ggplot(gsec_picmin_NB_merge,aes(x = diffs,y = contributing_prop)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~climate_var_F,scales = "free",strip.position = "top") +
  theme_bw() +
  theme() +
  labs(x = "Niche Breadth (Standardised by species mean)",
       y = "Proportional contributions to repeatability") +
  ggtitle("Niche Breadth - Standardised range relative to species mean")


# Final Figs --------------------------------------------------------------
pdf("figs/FigureSX_NB_and_GSEC_vs_repeatability.pdf",width=12,height = 10)
GSEC_fig + theme(axis.title = element_text(size = 16),
                 axis.text = element_text(size = 14),
                 axis.text.x = element_text(size = 14,angle = 45,hjust = 1))

NB1_fig + theme(axis.title = element_text(size = 16),
                axis.text = element_text(size = 14),
                axis.text.x = element_text(size = 14,angle = 45,hjust = 1),
                title = element_text(size = 18))

NB2_fig + theme(axis.title = element_text(size = 16),
                axis.text = element_text(size = 14),
                axis.text.x = element_text(size = 14,angle = 45,hjust = 1),
                title = element_text(size = 18))
dev.off()

# # Further analysis of difference in pvals by GSEC/Niche Breadth... --------
# # Take from the above and add in dataset, GSEC, NB
# repeatable_OG_pvals_merge = merge(repeatable_OG_pvals,gsec_picmin_NB_merge[,.(climate_var,dataset,species,GSEC,relative_breadth)],by = c("species","climate_var"),allow.cartesian=TRUE)
# 
# # For all species/climate, what's the average contributing pval...
# avg_pval_species_climate = repeatable_OG_pvals_merge[,.(avg_pval = mean(-log10(.SD$min_sdP_DS))),by = .(species,climate_var)]
# 
# # Merge with GSEC/NB
# avg_pval_species_climate = merge(avg_pval_species_climate,unique(gsec_picmin_NB_merge[,.(climate_var,species,dataset,GSEC,relative_breadth)]),
#                                  by = c("species","climate_var"))
# 
# # Plot each...
# ggplot(avg_pval_species_climate,aes(GSEC,avg_pval)) +
#   # facet_wrap(~species,scales = "free") + 
#   geom_point() +
#   geom_smooth(method = "lm") 
# 
# # # Plot each...
# # ggplot(avg_pval_species_climate,aes(relative_breadth,avg_pval)) +
# #   facet_wrap(~species,scales = "free") +
# #   geom_point() +
# #   geom_smooth(method = "lm") 



