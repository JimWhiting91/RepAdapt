# Script looks at association between WZA and recombination rate across subset of datasets...
lib = c("data.table","ggplot2","pbmcapply")
sapply(lib,library,character.only = T)
source("R/repadapt_functions.R")
source("R/PicMin.R")
n_cores = 6

# Which run name are we looking at?
run_name = "230321"
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"
pvals_file = paste0("outputs/GEA_res/run",run_name,"_",output_name,"_WZA_OG_pvals.rds")

# Fetch picmin results
picmin_outputs = readRDS(paste0("outputs/",output_name,"_picmin_results_doubleDS_pvals.rds"))
picmin_fdr = rbindlist(lapply(picmin_outputs,'[[',1))
OG_pvals = rbindlist(lapply(picmin_outputs,'[[',2))

# Fetch all the recombination specific pvals
rec_pval_files = grep(output_name,list.files('outputs/GEA_res/',pattern = "_RecRate_"),value = T)
rec_pvals = rbindlist(pbmclapply(rec_pval_files,function(x){
  # Get climate var
  clim_tmp = gsub('_pvals.rds','',gsub(paste0("run",run_name,"_",output_name,"_RecRate_"),'',x))
  # Read in data
  tmp = readRDS(paste0('outputs/GEA_res/',x)) |>
    mutate(climate_var = clim_tmp)
  tmp
},mc.cores = 4))

# Plot rec rate stuff -----------------------------------------------------

# Let's pick 4 random climate vars
set.seed(1000)
clim_random = sample(unique(rec_pvals$climate_var),4)

# Plot per gene pval vs rec rate pval
ggplot(rec_pvals[climate_var %in% clim_random,],aes(snp_count_pval,sd_epvalue)) + 
  ggdensity::geom_hdr()+
  facet_grid(dataset~climate_var) +
  theme_minimal()+
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 30))

pdf("figs/FigureSX_WZA_recombination_association.pdf",width = 12,height = 10)
ggplot(rec_pvals[climate_var %in% clim_random,],aes(rec_rate_pval,sd_epvalue)) + 
  ggdensity::geom_hdr()+
  facet_grid(dataset~climate_var) +
  theme_minimal()+
  theme(panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0),
        strip.text.x = element_text(angle = 30)) +
  labs(y = "WZA ep-value",x = "Recombination ep-value")
dev.off()


# Run PicMin over the recombination pvals ---------------------------------
# Remove duplicate species
rec_pvals_uniq = unique(rec_pvals[dataset %in% c("Arabidopsis_lyrata_Willi_PoolSeq",
                                                 "Arabidopsis_thaliana_Weigel_WGS_SWE",
                                                 "Arabis_alpina_Rellstab_PoolSeq",
                                                 "Helianthus_annuus_Todesco_Individual"),.(species,Orthogroup,rec_rate_pval)])

# Calc min P and DS correct
rec_minP = rec_pvals_uniq[,.(recrate_DS_pval = DunnSidak_pvals(rec_rate_pval)),by = .(species,Orthogroup)]

# Make sure we're testing the same orthogroups...
# and empiricise
rec_minP = rec_minP[Orthogroup %in% unique(picmin_fdr$Orthogroup),.(recrate_DS_epval = EmpiricalPs(recrate_DS_pval),
                                                                    Orthogroup),by = species]

# Run PicMin over these as we run it for GEA...
picmin_reps = 1e6
# Make our random tippet nulls
recrate_tippett_nulls = generateTippettNull(N = 4,numReps = picmin_reps)
# Random res for 4 species
# Make random data
random_input = GenerateRandomGenes(Ngenes = length(unique(picmin_fdr$Orthogroup)),
                                   Nsims = picmin_reps,
                                   Ngroups = 4)
# Run it through picmin
random_res = PicMin_bulk(pvals = random_input$pval,
                         groups = random_input$gene,
                         nulls = recrate_tippett_nulls)

# Now calculate observed recobination PicMin
to_test = table(rec_minP$Orthogroup)
to_test = names(to_test)[to_test == 4]
recrate_picmin_res = PicMin_bulk(pvals = rec_minP[Orthogroup %in% to_test,
                                                  recrate_DS_epval],
                                 groups = rec_minP[Orthogroup %in% to_test,
                                                   Orthogroup],
                                 nulls = recrate_tippett_nulls)

# Tidy up
recrate_picmin_res$picmin_p_adj = qvalue::empPvals(-recrate_picmin_res$picmin_p,-random_res$picmin_p)
recrate_picmin_res$picmin_fdr = p.adjust(recrate_picmin_res$picmin_p_adj,method = "fdr")
recrate_picmin_res$climate_var = 'Recombination'
colnames(recrate_picmin_res)[1] = 'Orthogroup'

# Any intersection here...
recrate_picmin_res[order(picmin_fdr),]
rec_picmin_hits = recrate_picmin_res[picmin_fdr < 0.5,Orthogroup]
climate_picmin_hits = unique(picmin_fdr[picmin_fdr < 0.5,Orthogroup])

# What are the tested props...
climate_prop = length(climate_picmin_hits) / length(unique(picmin_fdr$Orthogroup))
rec_prop = length(rec_picmin_hits) / length(unique(recrate_picmin_res$Orthogroup))
climate_prop * rec_prop * length(unique(recrate_picmin_res$Orthogroup))
sum(rec_picmin_hits %in% picmin_fdr[picmin_fdr < 0.5,Orthogroup])
