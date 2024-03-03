# Final PicMin Run
# Orthogroup-level convergence using PicMin
# This is a variant where the number of paralogs is corrected for by doing 2 empirical pvalue corrections with Dunn-Sidak
lib <- c("Rfast","qvalue","cowplot","poolr","mvmeta","qvalue","tidyr","ape","VGAM","ggExtra","pbmcapply","parallel","wCorr","data.table","ggplot2","viridis","ggridges","dplyr","readr")
sapply(lib,library,character.only=T)

# Function Library --------------------------------------------------------
source("R/PicMin.R")
DunnSidak_pvals = function(pval_vector){
  1 - (1 - min(pval_vector))^length(pval_vector)
}

stouffer_pvals = function(pval_vector){
  sum(qnorm(1 - pval_vector)) / sqrt(length(pval_vector))
}

stouffer_pvals_noMin = function(pval_vector){
  pval_vector2 = pval_vector[-which.min(pval_vector)]
  sum(qnorm(1 - pval_vector2)) / sqrt(length(pval_vector2))
}

############################################################################
# What GEA results are we looking at?
run_name = "230321"
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"
pvals_file = paste0("outputs/GEA_res/run",run_name,"_",output_name,"_WZA_OG_PerGene_pvals.rds")
seedN = 1000 # Seed for reproducible runs
picmin_reps = 1e6 # The number of permutations used to test significance of picmin pvals
n_cores = 6 # Number of cpu
orthogroup_cutoff <- 20 # Minimum number of species in an OG
maxParalogN = 10 # What's the maximum number of paralogs in a given orthogroup, this is already filtered but needed to set up some inputs

# Fetch the processed Orthogroup Pvals and only retain the pvalues from one climate variable
OG_pergene_pvals = readRDS(pvals_file)

# Define the brassica species, and set a new output file name
output_name2 = 'brassica_ONLY'
brassica_species = grep('Arab|Cardam|Boec|Capsell',unique(OG_pergene_pvals$mean_temp$species),
                        value = T)
orthogroup_cutoff = length(brassica_species)
# And filter each of the pval files
for(i in 1:length(OG_pergene_pvals)){
  OG_pergene_pvals[[i]] = OG_pergene_pvals[[i]][species %in% brassica_species,]
}

# Use the information from the first set of pvals to define some variables
dummy_data = unique(OG_pergene_pvals[[1]][,.(species,Orthogroup)])
OG_counts = table(dummy_data$Orthogroup)
test_OG = names(OG_counts)[OG_counts >= orthogroup_cutoff]


# Get PicMin nulls ready --------------------------------------------------
# Set up the Tippett Nulls
tippett_nulls = pbmclapply(orthogroup_cutoff:length(unique(dummy_data$species)),generateTippettNull,geneN = 10000,seedN = 1000,numReps = 1000000,mc.cores = n_cores)
names(tippett_nulls) = as.character(orthogroup_cutoff:length(unique(dummy_data$species)))

# Set up empirical null under randomness
picmin_nulls = lapply(orthogroup_cutoff:length(unique(dummy_data$species)),function(x){
  message(paste0(">>> Making random picmin results for ",x," species..."))
  # Make random data
  random_input = GenerateRandomGenes(Ngenes = length(unique(dummy_data$Orthogroup)),
                                     Nsims = picmin_reps,
                                     Ngroups = x)
  # Run it through picmin
  random_res = PicMin_bulk(pvals = random_input$pval,
                           groups = random_input$gene,
                           nulls = tippett_nulls[[as.character(x)]])
  random_res$picmin_p
})
names(picmin_nulls) = as.character(orthogroup_cutoff:length(unique(dummy_data$species)))

# Plot the random nulls for demonstration purposes...
plot_nulls = data.table(pvals = unlist(picmin_nulls),
                        species = rep(names(picmin_nulls),each = picmin_reps))

# pdf("figs/FigureSX_random_picmin_pval_distributions.pdf",width = 8,height = 6)
ggplot(plot_nulls,aes(pvals)) +
  geom_histogram(breaks = seq(0, 1, by = 1/100)) +
  facet_wrap(~species) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  labs(y = "Count",x = "PicMin p-value")
# dev.off()

# Run PicMin in Bulk ------------------------------------------------------
# Run over each climate separately
all_picmin_res = lapply(names(OG_pergene_pvals),function(focal_climate){
  message(paste0("Running ",focal_climate))
  
  # Filter Orthogroups for those we're testing...
  OG_pvals_climate = OG_pergene_pvals[[focal_climate]]
  
  # Take empirical pval
  OG_pvals_climate = OG_pvals_climate[,.(Orthogroup,
                                         combined_sd_pvalue,
                                         epval = EmpiricalPs(combined_sd_pvalue)),by = species]
  
  OG_counts = table(unique(OG_pvals_climate[,.(species,Orthogroup)])$Orthogroup)
  OG_pvals_test = OG_pvals_climate[Orthogroup %in% names(OG_counts)[OG_counts >= orthogroup_cutoff],]
  
  # Run PicMin
  message(">>> Adjusting p-values")
  OG_pvals_test_DS = OG_pvals_test[,.(epval_DS = DunnSidak_pvals(epval),
                                      Ngenes_per_species = nrow(.SD)),by = .(species,Orthogroup)]
  OG_pvals_test_DS = OG_pvals_test_DS[,.(Orthogroup,epval_DS,Ngenes_per_species,
                                         epval_final = EmpiricalPs(epval_DS)),by = species]
  # Run PicMin
  message(">>> Running PicMin and doing corrections")
  picmin_res = PicMin_bulk(pvals = OG_pvals_test_DS$epval_final,
                           groups = OG_pvals_test_DS$Orthogroup,
                           nulls = tippett_nulls$`7`)
  picmin_res$climate_var = focal_climate
  
  # Adjust for the shape of the null distribution
  picmin_res = rbindlist(lapply(unique(picmin_res$testN),function(x){
    tmp = picmin_res[testN == x,]
    tmp$picmin_p_adj = qvalue::empPvals(-tmp$picmin_p,-picmin_nulls[[as.character(x)]])
    tmp
  }))
  
  # FDR correct and tidy
  picmin_res$picmin_fdr = p.adjust(picmin_res$picmin_p_adj,method = 'fdr')
  colnames(picmin_res)[1] = 'Orthogroup'
  
  # Run Stouffers as a comparison
  stouffer_res = OG_pvals_test_DS[,.(stoufferZ = stouffer_pvals(epval_final)),by = Orthogroup]
  stouffer_res$stouffer_pval = pnorm(stouffer_res$stoufferZ,lower.tail = F)
  stouffer_res$stouffer_fdr = p.adjust(stouffer_res$stouffer_pval,'fdr')
  
  # Return merged results
  OG_pvals_test_DS$climate_var = focal_climate
  return(list(picmin_res = merge(picmin_res,stouffer_res),
              tested_pvals = OG_pvals_test_DS))
})

# Save these and export
saveRDS(all_picmin_res,
        paste0("outputs/",output_name2,"_picmin_results_doubleDS_pvals.rds"))


# Compare these to original RAOs ------------------------------------------
# Fetch the original picmin res
picmin_outputs = readRDS(paste0("outputs/",output_name,"_picmin_results_doubleDS_pvals.rds"))
full_picmin_res = rbindlist(lapply(picmin_outputs,'[[',1))
# what OG were tested
original_tested_OG = unique(full_picmin_res$Orthogroup)

# Also fetch the Athal map
OG_map_Athal = data.table(readRDS(paste0("data/OG_map_Athal_",output_name,".rds")))

# Fetch the brassica signif...
brassica_signif = lapply(all_picmin_res,'[[',1) |>
  rbindlist()
brassica_signif[picmin_fdr < 0.5,][order(picmin_fdr),]


# Add a variable for whether or not it was previously tested...
brassica_signif$previously_tested = ifelse(brassica_signif$Orthogroup %in% original_tested_OG,
                                           'Full Dataset','Brassica Dataset')
# And another variable for whether they are UNIQUE to brassica only...
all_pvals = readRDS(pvals_file)
# Which orthogroups ARE IN SPECIES THAT ARE NOT IN BRASSICA
noBrassica_OG = rbindlist(all_pvals)[!species %in% brassica_species,Orthogroup] |>
  unique()
# Which orthogroups ARE IN SPECIES THAT ARE NOT IN BRASSICA
Brassica_OG = rbindlist(all_pvals)[species %in% brassica_species,Orthogroup] |>
  unique()
# Which Orthogroups are then UNIQUE to BRASSICA
brassica_uniq_OG = Brassica_OG[!Brassica_OG %in% noBrassica_OG]
# Add the variable
brassica_signif$brassica_unique = ifelse(brassica_signif$Orthogroup %in% brassica_uniq_OG,
                                         'Brassica unique','Not unique')

# General OG counts
length(unique(brassica_signif$Orthogroup))
length(unique(brassica_signif$Orthogroup)[unique(brassica_signif$Orthogroup) %in% original_tested_OG])
length(unique(brassica_signif$Orthogroup)[unique(brassica_signif$Orthogroup) %in% brassica_uniq_OG])

# what proportions are these?
table(brassica_signif$previously_tested) / nrow(brassica_signif)
table(brassica_signif[picmin_fdr < 0.5,previously_tested]) / nrow(brassica_signif[picmin_fdr < 0.5,])

table(brassica_signif$brassica_unique) / nrow(brassica_signif)
table(brassica_signif[picmin_fdr < 0.5,brassica_unique]) / nrow(brassica_signif[picmin_fdr < 0.5,])

# Plot them side by side
brassica_picmin_test = cowplot::plot_grid(
  ggplot(brassica_signif,aes(fill = previously_tested,x = picmin_p_adj)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c('Full Dataset' = 'orange3',
                                 'Brassica Dataset' = 'blue4'))  +
    labs(y = 'Density',
         x = 'PicMin P-Value',
         fill = 'Orthogroup Status') +
    theme_classic() +
    theme(legend.position = 'top',
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12)),
  ggplot(brassica_signif,aes(fill = brassica_unique,x = picmin_p_adj)) +
    geom_density(alpha = 0.5) +
    scale_fill_manual(values = c('Not unique' = 'orange3',
                                 'Brassica unique' = 'blue4')) +
    labs(y = 'Density',
         x = 'PicMin P-Value',
         fill = 'Orthogroup Status') +
    theme_classic() +
    theme(legend.position = 'top',
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12)),
  nrow = 1,align = 'h',axis = 'tblr'
)
brassica_picmin_test

# Save these results
pdf('figs/FigureSX_brassica_only_picmin_test.pdf',width = 12,height = 4)
brassica_picmin_test
dev.off()

# Also what about brassica unique vs not tested in the main analysis but not unique
ggplot(brassica_signif[previously_tested != 'Full Dataset',],aes(fill = brassica_unique,x = picmin_p_adj)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c('Not unique' = 'orange3',
                               'Brassica unique' = 'blue4')) +
  labs(y = 'Density',
       x = 'PicMin P-Value',
       fill = 'Orthogroup Status') +
  theme_classic() +
  theme(legend.position = 'top',
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12))

# Just also test quickly for differences in gene duplication between the groups...
OG_paralogs = OG_pergene_pvals[[1]][,.(paralog_N = nrow(.SD)),by = .(species,Orthogroup)]
OG_paralogs = OG_paralogs[,.(avg_paralog_N = mean(paralog_N)),by = Orthogroup]

# Attach these and estimate diffs
merge(brassica_signif,OG_paralogs,all.x = T) |>
subset(climate_var == 'mean_temp') |>
  ggplot(aes(x = previously_tested,y = avg_paralog_N)) +
  geom_boxplot() +
  scale_y_continuous(trans = 'log')

# brassica_signif$climate_type = ifelse(grepl('prec',brassica_signif$climate_var),
#                                       'Precipitation','Temperature')
# ggplot(brassica_signif,aes(fill = climate_type,x = picmin_p_adj)) +
#   geom_density(alpha = 0.5) +
#   facet_wrap(~previously_tested)

