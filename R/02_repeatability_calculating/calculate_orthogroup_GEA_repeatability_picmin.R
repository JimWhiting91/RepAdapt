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

# Use the information from the first set of pvals to define some variables
dummy_data = unique(OG_pergene_pvals[[1]][,.(species,Orthogroup)])
OG_counts = table(dummy_data$Orthogroup)
test_OG = names(OG_counts)[OG_counts >= orthogroup_cutoff]


# Plot distributions of tested vs untested pvals --------------------------
all_pvals = rbindlist(OG_pergene_pvals)
all_pvals = all_pvals[,.(Orthogroup,
                         epval = EmpiricalPs(combined_sd_pvalue)),by = .(species,climate_var)]
all_pvals$tested = "Tested"
all_pvals[!Orthogroup %in% test_OG,tested := "Not Tested"]

# Plot
pdf("figs/FigureSX_tested_nottested_pvalue_distributions.pdf",width = 12,height = 12)
ggplot(all_pvals,aes(epval,fill = tested)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~species) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12,angle = 45,hjust = 1),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 12)) +
  labs(x = "GEA ep-value",y = "Density",fill = "Tested for\nRepeatability?")
dev.off()


# Get PicMin nulls ready --------------------------------------------------
# Set up the Tippett Nulls
tippett_nulls = pbmclapply((orthogroup_cutoff - 1):length(unique(dummy_data$species)),generateTippettNull,geneN = 10000,seedN = 1000,numReps = 1000000,mc.cores = n_cores)
names(tippett_nulls) = as.character((orthogroup_cutoff - 1):length(unique(dummy_data$species)))

# Set up empirical null under randomness
picmin_nulls = lapply((orthogroup_cutoff - 1):length(unique(dummy_data$species)),function(x){
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
names(picmin_nulls) = as.character((orthogroup_cutoff - 1):length(unique(dummy_data$species)))

# Plot the random nulls for demonstration purposes...
plot_nulls = data.table(pvals = unlist(picmin_nulls),
                        species = rep(names(picmin_nulls),each = picmin_reps))

pdf("figs/FigureSX_random_picmin_pval_distributions.pdf",width = 8,height = 6)
ggplot(plot_nulls,aes(pvals)) +
  geom_histogram(breaks = seq(0, 1, by = 1/100)) +
  facet_wrap(~species) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  labs(y = "Count",x = "PicMin p-value")
dev.off()

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
  # There is a catch here for isothermality because there is no GEA data from one dataset
  if(focal_climate == 'isothermality'){
    OG_pvals_test = OG_pvals_climate[Orthogroup %in% names(OG_counts)[OG_counts >= (orthogroup_cutoff - 1)],]
  } else {
    OG_pvals_test = OG_pvals_climate[Orthogroup %in% names(OG_counts)[OG_counts >= orthogroup_cutoff],]
  }
  
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
                           nulls = tippett_nulls)
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
        paste0("outputs/",output_name,"_picmin_results_doubleDS_pvals.rds"))


# Repeat PicMin, but randomise pvals prior to paralog DS ------------------
randomN = 1000

# Shuffle all the orthogroups randomN times
# Filter Orthogroups for those we're testing...
OG_pvals_all = rbindlist(OG_pergene_pvals)
OG_counts = table(unique(OG_pvals_all[,.(species,Orthogroup)])$Orthogroup)
to_test = names(OG_counts)[OG_counts >= orthogroup_cutoff]

all_random_picmin_res = rbindlist(lapply(1:randomN,function(x){
  set.seed(x)
  
  # Take empirical pval
  OG_pvals_climate = OG_pvals_all[Orthogroup %in% to_test,
                                  .(Orthogroup,random_epval = sample(EmpiricalPs(combined_sd_pvalue))),
                                  by = .(species,climate_var)]
  
  # Do second DS
  message(">>> Adjusting p-values")
  OG_pvals_test = OG_pvals_climate[Orthogroup %in% to_test,]
  OG_pvals_test_DS = OG_pvals_test[,.(random_epval_DS = DunnSidak_pvals(random_epval),
                                      Ngenes_per_species = nrow(.SD)),by = .(species,Orthogroup,climate_var)]
  OG_pvals_test_DS = OG_pvals_test_DS[,.(Orthogroup,random_epval_DS,Ngenes_per_species,
                                         random_epval_final = EmpiricalPs(random_epval_DS)),by = .(species,climate_var)]
  
  OG_pvals_test_DS$random_iter = x
  
  # Set up test groups
  OG_pvals_test_DS$test_group = paste0(OG_pvals_test_DS$Orthogroup,":",OG_pvals_test_DS$climate_var)
  
  # Run PicMin
  random_picmin_res = PicMin_bulk(OG_pvals_test_DS$random_epval_final,
                                  OG_pvals_test_DS$test_group,
                                  tippett_nulls)
  
  # Adjust for the shape of the null distribution
  random_picmin_res = rbindlist(lapply(unique(random_picmin_res$testN),function(x){
    tmp = random_picmin_res[testN == x,]
    tmp$picmin_p_adj = qvalue::empPvals(-tmp$picmin_p,-picmin_nulls[[as.character(x)]])
    tmp
  }))
  
  random_picmin_res$iteration = x
  random_picmin_res
}))


# Also save these...
saveRDS(all_random_picmin_res,
        paste0("outputs/",output_name,"_randomised_pvals_picmin_res.rds"))

# all_random_picmin_res = readRDS(
#   paste0("outputs/",output_name,"_randomised_pvals_picmin_res.rds")
# )

# How many RAO with FDR < 0.5 under random --------------------------------
# Fetch back the results
all_random_picmin_res = readRDS(
  paste0("outputs/",output_name,"_randomised_pvals_picmin_res.rds")
)

# Split by variable
# FDR-adjust within each iteration and within each clim var...
clim_specific_random = lapply(names(OG_pergene_pvals),function(clim_var){
  print(paste0('>>> Starting ',clim_var))
  
  clim_tmp = all_random_picmin_res[grepl(clim_var,test_group),] |>
    subset(!grepl(paste0(clim_var,'_'),test_group)) |>
    tidyr::separate(test_group,sep = ':',into = c('Orthogroup','climate_var')) |>
    as.data.table()
  
  # Filter duplicate temp names
  clim_tmp = clim_tmp[climate_var == clim_var,]
  print(paste0('Number of rows: ',nrow(clim_tmp)))
  
  # Calc FDR in climate vars and iterations
  random_FDR = clim_tmp[,.(picmin_fdr = p.adjust(.SD$picmin_p_adj,method = 'fdr'),
                           Orthogroup),
                        by = .(iteration)]
  # How many with FDR < 0.5 to be expected
  exp_fdr = random_FDR[,.(exp_05 = sum(picmin_fdr < 0.5),
                          exp_04 = sum(picmin_fdr < 0.4),
                          exp_03 = sum(picmin_fdr < 0.3),
                          exp_02 = sum(picmin_fdr < 0.2),
                          exp_01 = sum(picmin_fdr < 0.1)),by = iteration]
  # We also need to make a list per iteration of all significant orthogroups for OG collapsing...
  signif_OG = random_FDR[picmin_fdr < 0.5,]
  exp_fdr$climate_var = clim_var
  signif_OG$climater_var = clim_var
  return(list(exp_fdr,signif_OG))
})

# Save these
saveRDS(clim_specific_random,
        paste0("outputs/",output_name,"_randomised_pvals_picmin_res_final_FDRs.rds"))

# # Split out into raw numbers and OG
# random_FDR_counts = lapply(clim_specific_random,'[[',1) |>
#   rbindlist()
# all_iterations = random_FDR_counts[,.(sum_fdr_05 = sum(exp_05),
#                                       sum_fdr_04 = sum(exp_04),
#                                       sum_fdr_03 = sum(exp_03),
#                                       sum_fdr_02 = sum(exp_02),
#                                       sum_fdr_01 = sum(exp_01)),by = iteration]
# 
# 
# hist(all_iterations$sum_fdr_05)
# hist(all_iterations$sum_fdr_04)
# hist(all_iterations$sum_fdr_03)
# hist(all_iterations$sum_fdr_02)
# hist(all_iterations$sum_fdr_01)
# 
# 
# mean(all_iterations$sum_fdr_05)
# mean(all_iterations$sum_fdr_04)
# mean(all_iterations$sum_fdr_03)
# mean(all_iterations$sum_fdr_02)
# mean(all_iterations$sum_fdr_01)

