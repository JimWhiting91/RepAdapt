# This script takes the results from GEA scans, and condenses them all down into a single table that shows for a given OG, Dataset, Climate Var, a single pval
lib <- c("poolr","regioneR","Rfast","mvmeta","qvalue","tidyr","ggtree","ape","VGAM","ggExtra","WeMix","pbmcapply","parallel","wCorr","data.table","ggplot2","viridis","ggridges","dplyr","patchwork","readr")
sapply(lib,library,character.only=T)

# What GEA results are we looking at and processing...
run_name = "230321"
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"

# Get the map of OG - GEA genes...
OG_map = readRDS(paste0("data/OF_OG_gea_gene_map_",output_name,".rds"))

# And the picmin results
picmin_outputs = readRDS(paste0("outputs/",output_name,"_picmin_results_doubleDS_pvals.rds"))
picmin_fdr = rbindlist(lapply(picmin_outputs,'[[',1))

# From the map, what proportion of genes are actually tested
OG_map$tested = "Not Tested"
OG_map[Orthogroup %in% unique(picmin_fdr$Orthogroup), tested := "Tested"]

# Get gene props...
tested_prop = na.omit(OG_map)[,.(test_prop = nrow(.SD[tested == "Tested",])/nrow(.SD)),by = .(genome)][order(test_prop),]

# What are single copy rates?
paralog_counts = na.omit(OG_map)[,data.table(table(table(.SD$Orthogroup))),by = genome]
single_copy = paralog_counts[V1 == 1,][order(genome),]
total = paralog_counts[,.(total = sum(N)),by = genome][order(genome)]
single_copy = merge(single_copy,total)                     
sort(single_copy$N / single_copy$total)

# Make a supp-table...
out = merge(tested_prop,single_copy[,.(N,total,genome)],by = "genome")
colnames(out) = c("Genome","Tested Proportion","Single Copy N","Total N")
out$`Single Copy Proportion` = round(out$`Single Copy N`/out$`Total N`,3)
out$`Tested Proportion` = round(out$`Tested Proportion`,3)

# Also add the data from the comparative stats from orthofinder...
percent_assigned = data.table(Genome = sort(unique(out$Genome)),
                              `% Genes Assigned to Orthogroups` = c(91.2,96.0,99.4,96.7,91.7,96.9,99.4,96.7,86.7,88.0,85.7,95.0,92.9,93.0,91.9,89.4,96.7,96.0))

# Save this
write.csv(merge(out,percent_assigned),"tables/TableSX_orthogroup_summary_stats.csv",
          row.names = F,quote = F)
