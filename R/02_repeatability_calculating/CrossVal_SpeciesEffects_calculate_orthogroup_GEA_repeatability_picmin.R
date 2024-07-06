# PicMin Repeatability Tests
# Orthogroup-level convergence using PicMin
# This is a variant of the PicMin process that does cross-validation to explore effects of individual datasets
# PicMin is run for all tested orthogroups, but each species is excluded once from all tests
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


# Get PicMin nulls ready --------------------------------------------------
# Set up the Tippett Nulls
# We nulls for -1 sets given we will be excluding a species from CV
tippett_nulls = pbmclapply((orthogroup_cutoff - 2):(length(unique(dummy_data$species)) - 1),generateTippettNull,geneN = 10000,seedN = 1000,numReps = 1000000,mc.cores = n_cores)
names(tippett_nulls) = as.character((orthogroup_cutoff - 2):(length(unique(dummy_data$species)) - 1))

# Set up empirical null under randomness
picmin_nulls = lapply((orthogroup_cutoff - 2):(length(unique(dummy_data$species)) - 1),function(x){
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
names(picmin_nulls) = as.character((orthogroup_cutoff - 2):(length(unique(dummy_data$species)) - 1))

# Plot the random nulls for demonstration purposes...
plot_nulls = data.table(pvals = unlist(picmin_nulls),
                        species = rep(names(picmin_nulls),each = picmin_reps))

ggplot(plot_nulls,aes(pvals)) +
  geom_histogram(breaks = seq(0, 1, by = 1/100)) +
  facet_wrap(~species) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14)) +
  labs(y = "Count",x = "PicMin p-value")

# Run PicMin in Bulk ------------------------------------------------------
# For the cross-validations, we need to exclude one species from the entire analysis
# and re-run
true_picmin_outputs = readRDS(paste0("outputs/",output_name,"_picmin_results_doubleDS_pvals.rds"))
true_picmin_fdr = rbindlist(lapply(picmin_outputs,'[[',1))

# To make sure we are testing the same orthogroups, read back in the real results and filter
# for only those tested orthogroups...
species_cross_vals = lapply(unique(unique(OG_pergene_pvals[[1]]$species)),function(spec){
  message(paste0('>>> Excluding ',spec))
  
  # Run over each climate separately
  all_picmin_res = lapply(names(OG_pergene_pvals),function(focal_climate){
    message(paste0("Running ",focal_climate))
    
    # Filter Orthogroups for those we're testing...
    OG_pvals_climate = OG_pergene_pvals[[focal_climate]]
    
    # Which OG were tested for this variable??
    true_tested_OG = true_picmin_fdr[climate_var == focal_climate,Orthogroup]
    # EXCLUDE THE FOCAL SPECIES AND INCLUDE ONLY TRUE TESTED OG
    OG_pvals_climate = OG_pvals_climate[species != spec & Orthogroup %in% true_tested_OG,]
    
    # Take empirical pval
    OG_pvals_climate = OG_pvals_climate[,.(Orthogroup,
                                           combined_sd_pvalue,
                                           epval = EmpiricalPs(combined_sd_pvalue)),by = species]
    
    # OG_counts = table(unique(OG_pvals_climate[,.(species,Orthogroup)])$Orthogroup)
    # # There is a catch here for isothermality because there is no GEA data from one dataset
    # # For cross-vals these are one less than usual
    # if(focal_climate == 'isothermality'){
    #   OG_pvals_test = OG_pvals_climate[Orthogroup %in% names(OG_counts)[OG_counts >= (orthogroup_cutoff - 2)],]
    # } else {
    #   OG_pvals_test = OG_pvals_climate[Orthogroup %in% names(OG_counts)[OG_counts >= orthogroup_cutoff - 1],]
    # }
    OG_pvals_test = OG_pvals_climate
    
    # Run PicMin
    # message(">>> Adjusting p-values")
    OG_pvals_test_DS = OG_pvals_test[,.(epval_DS = DunnSidak_pvals(epval),
                                        Ngenes_per_species = nrow(.SD)),by = .(species,Orthogroup)]
    OG_pvals_test_DS = OG_pvals_test_DS[,.(Orthogroup,epval_DS,Ngenes_per_species,
                                           epval_final = EmpiricalPs(epval_DS)),by = species]
    # Run PicMin
    # message(">>> Running PicMin and doing corrections")
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
    return(merge(picmin_res,stouffer_res))
  }) |>
    rbindlist()
  
  # Which species was dropped?
  all_picmin_res$focal_species = spec
  return(all_picmin_res)
})

# Save these and export
saveRDS(species_cross_vals,
        paste0("outputs/",output_name,"_picmin_results_SpeciesCrossValidation.rds"))


