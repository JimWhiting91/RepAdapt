#### Functional analysis of genes identified as convergent via PicMin analysis
lib = c("cowplot","ggridges","viridis","tidyr","ape","DHARMa","sjPlot","qvalue","ggplot2","pbmcapply","biomaRt","org.At.tair.db","ExpressionAtlas","S4Vectors", "IRanges", "GenomicRanges", "SummarizedExperiment","readr","data.table")
sapply(lib,library,character.only=T)
n_cores = 6

species_genome_map = data.table(species = c("Pseudotsuga menziesii",
                                            "Picea abies","Picea obovata","Picea glaucaxengelmannii",
                                            "Pinus sylvestris","Pinus contorta",
                                            "Panicum hallii",
                                            "Amaranthus tuberculatus",
                                            "Helianthus annuus","Helianthus argophyllus","Helianthus petiolaris",
                                            "Medicago truncatula",
                                            "Eucalyptus albens","Eucalyptus sideroxylon","Eucalyptus magnificata",
                                            "Quercus petraea",
                                            "Populus deltoides",
                                            "Populus trichocarpa",
                                            "Populus tremula",
                                            "Arabis alpina",
                                            "Capsella rubella","Capsella bursa-pastoris",
                                            "Boechera stricta",
                                            "Arabidopsis thaliana",
                                            "Arabidopsis halleri","Cardamine resedifolia",
                                            "Arabidopsis lyrata"),
                                genome = c("Pmenziesii",
                                           "Pabies","Pabies","Pabies",
                                           "Ptaeda","Ptaeda",
                                           "Phallii",
                                           "Atuberculatus",
                                           "Hannuus","Hannuus","Hannuus",
                                           "Mtruncatula",
                                           "Egrandis","Egrandis","Egrandis",
                                           "Qpetraea",
                                           "Pdeltoides",
                                           "Ptrichocarpa",
                                           "Ptremula",
                                           "Aalpina",
                                           "Crubella","Crubella",
                                           "Bstricta",
                                           "Athaliana",
                                           "Ahalleri","Ahalleri",
                                           "Alyrata"))

# Functions ---------------------------------------------------------------
calcTau <- function(input_vec){
  input_vec2 = as.numeric(input_vec) * 10
  # input_vec2[input_vec2 < 1] <- 0
  TPM_max = max(input_vec2)
  if(TPM_max >= 1){
    sum(1 -  input_vec2/TPM_max) / (length(input_vec2) - 1)
  } else {
    return(as.numeric(NA))
  }
}
calcLogTau <- function(input_vec){
  input_vec2 = as.numeric(input_vec) * 10
  # input_vec2[input_vec2 < 1] <- 0
  TPM_max = max(input_vec2)
  if(TPM_max >= 1){
    norm_vec = log(input_vec2)/log(TPM_max)
    norm_vec[!is.finite(norm_vec)] = 0
    return(sum(1 - norm_vec) / (length(input_vec) - 1))
  } else {
    return(as.numeric(NA))
  }
}

condenseToOGPvals <- function(data,OG_col,pval_col){
  
  data$Orthogroup = data.frame(data)[,OG_col]
  data$tmp_p = data.frame(data)[,pval_col]
  
  # Group summarise
  OG_minP = data[,.(minP = min(.SD$tmp_p),
                    N = nrow(.SD)),by = Orthogroup]
  
  # Convert to DS corrected
  OG_minP$minP_DS =  1 - (1 - OG_minP$minP)^OG_minP$N
  
  # Add epvalue
  OG_minP$minP_DS_epvalue = empPvals(-OG_minP$minP_DS,-OG_minP$minP_DS)
  OG_minP
}

stoufferZ <- function(Z_vector){
  sum(Z_vector) / sqrt(length(Z_vector))
}


# Which run name are we looking at?
run_name = "230321"
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"

# Fetch picmin results
picmin_outputs = readRDS(paste0("outputs/",output_name,"_picmin_results_doubleDS_pvals.rds"))
picmin_fdr = rbindlist(lapply(picmin_outputs,'[[',1))
OG_pvals = rbindlist(lapply(picmin_outputs,'[[',2))

# Split the picmin results into per orthogroup strongest repeatability deciles
picmin_deciles = picmin_fdr[,.(PicMin_minP = min(picmin_p_adj)),by = Orthogroup]
# And exclude any orthogroups that were only tested once...
OG_counts = table(unique(OG_pvals[,.(Orthogroup,climate_var)])$Orthogroup)
picmin_deciles = picmin_deciles[Orthogroup %in% names(OG_counts)[OG_counts == max(OG_counts)]]
picmin_deciles$picmin_decile = cut_number(picmin_deciles$PicMin_minP, n = 10)
  
# Fetch OG - Athal gene map
OG_map_Athal = data.table(readRDS(paste0("data/OG_map_Athal_",output_name,".rds")))
# Fetch OG - Mtrunc gene map
OG_map_Mtrunc = data.table(readRDS(paste0("data/OG_map_Mtrunc_",output_name,".rds")))

# Prepare Expression Atlas Tau Scores -------------------------------------
# Also look at a more tissue specific one...
tissue_expression_atlas = data.table(read_tsv("data/E-MTAB-7978-query-results_all_tissue.tpms.tsv",comment = "#",skip = 4))

# Filter tissues for adult only, and identify unique tissues
all_tissues = grep("Gene",colnames(tissue_expression_atlas),invert=T,value=T)
tissue_types = data.table(variable = all_tissues,
                          tissue = sapply(strsplit(all_tissues,", "),'[[',2))
colnames(tissue_expression_atlas)[1] <- "Gene_ID"

# Transform the expression data to long-form...
tissue_expression_long = data.table(melt(tissue_expression_atlas,id.vars = "Gene_ID",measure.vars = grep("Gene",colnames(tissue_expression_atlas),invert = T,value = T))) 
tissue_expression_long = merge(tissue_expression_long,tissue_types)
tissue_expression_long$value = as.numeric(tissue_expression_long$value)

# Replace NA with 0...
tissue_expression_long$value[is.na(tissue_expression_long$value)] <- 0

# Summarise this within tissue types...
tissue_expression_sum = tissue_expression_long[,.(tpm_mean = mean(value),
                                                  tpm_max = max(value)),by = .(Gene_ID,tissue)]

# Now estimate tau using our tau10 estimate...
tissue_expression_tau = tissue_expression_sum[,.(tau_mean = calcLogTau(tpm_mean),
                                                 tau_max = calcLogTau(tpm_max),
                                                 tau_nolog = calcTau(tpm_mean),
                                                 max_tpm = max(tpm_mean)),
                                              by = Gene_ID]

# Fix stupid 1s...
tissue_expression_tau[tau_mean > 1,tau_mean := 1]
tissue_expression_tau[tau_max > 1,tau_max := 1]

# tissue_expression_tau = tissue_expression_sum[,.(tau_mean = sum(.SD$tpm_max)),by = Gene_ID]
tissue_expression_tau$TAIR_gene = tissue_expression_tau$Gene_ID

# # Save this
# saveRDS(tissue_expression_tau,
#         "outputs/Athal_tau_data.rds")


# Merge with the Orthogroups from OG map
pergene_tau_merge = merge(data.table(OG_map_Athal)[,.(TAIR_gene,Orthogroup)],tissue_expression_tau,by="TAIR_gene")
# Filter for Orthogroups that we actually test for convergence...
pergene_tau_merge = pergene_tau_merge[Orthogroup %in% unique(picmin_fdr$Orthogroup),]

# Convert tau scores into empirical pvals then Z scores, correcting through DS for Orthogroup...
# Looking for low tau...
pergene_tau_merge$tau_p1 = empPvals(-pergene_tau_merge$tau_nolog,-pergene_tau_merge$tau_nolog)
OG_tau1 = condenseToOGPvals(pergene_tau_merge,"Orthogroup","tau_p1")
OG_tau1$tau_Z = qnorm(OG_tau1$minP_DS_epvalue,lower.tail = F)


# We consider convergent OG as those with FDR < fdr_cutoff
# fdr_cutoff = 0.5
# convergent_OG = unique(picmin_fdr[picmin_fdr < fdr_cutoff,Orthogroup])
convergent_OG = unique(picmin_fdr[picmin_p_adj < 0.005,Orthogroup])

# Find the tau scores associated with convergent OG
hist(OG_tau1[Orthogroup %in% convergent_OG,tau_Z])
mean(OG_tau1[Orthogroup %in% convergent_OG,tau_Z])

# Merge with picmin deciles and estimate stouffers across
OG_tau1 = merge(OG_tau1,picmin_deciles)
picmin_decile_Zscores = OG_tau1[,.(stouffersZ = stoufferZ(tau_Z)),by = picmin_decile][order(picmin_decile),]
picmin_decile_Zscores$variable = 'Tau (Specificity)'

# We can also look at how specificity varies for repeatability of specific climate variables --------
# Run over climate vars
climate_vars = unique(picmin_fdr$climate_var)
tau_climate_stouffers = rbindlist(lapply(climate_vars,function(clim){
  
  convergent_OG = unique(picmin_fdr[picmin_p_adj < 0.005 & climate_var == clim,Orthogroup])
  
  # Find the tau scores associated with convergent OG
  clim_Z = OG_tau1[Orthogroup %in% convergent_OG,tau_Z]
  clim_stouffer = sum(clim_Z)/sqrt(length(clim_Z))
  
  # Output...
  out = data.table(climate_var = clim,
                   tested_N = length(convergent_OG),
                   tau_stouffer = clim_stouffer,
                   pval = pnorm(clim_stouffer,lower.tail = F))
}))
tau_climate_stouffers[order(-tau_stouffer),]


# Is this intermediate pleiotropy? ----------------------------------------
# Split into deciles and check which deciles are 'enriched'
OG_tau1$convergent = "No"
OG_tau1[Orthogroup %in% unique(picmin_fdr[picmin_p_adj < 0.005,Orthogroup]),convergent := "Yes"]

# What are the deciles...
OG_tau1$deciles = cut_number(OG_tau1$tau_Z,10)

# Plot these...
tau_deciles = data.table(obs = table(OG_tau1$deciles) / nrow(OG_tau1),
                         convergent = table(OG_tau1[convergent == "Yes",deciles]) / nrow(OG_tau1[convergent == "Yes",]))
tau_deciles$decile = factor(1:10,levels = 1:10)

tau_decile_plot = melt(tau_deciles[,.(decile,obs.N,convergent.N)]) |>
  ggplot(aes(y = value, x = factor(decile),fill = variable)) +
  geom_bar(stat = 'identity',position = 'dodge') +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        title = element_text(size = 16)) +
  labs(y = "Proportion",x = "Decile",fill = "Repeated?") +
  scale_fill_discrete(breaks = c("obs.N","convergent.N"),
                      labels = c("Expected","Repeated")) +
  ggtitle("Tissue specificity - Tau")
tau_decile_plot

# Plot examples of increasing specificity ---------------------------------
specificity_examples = sapply(quantile(pergene_tau_merge$tau_nolog,seq(0.1,0.9,0.2),na.rm = T),function(x){
  diffs = abs(pergene_tau_merge$tau_nolog - x)
  pergene_tau_merge[which(diffs == min(diffs,na.rm = T)),TAIR_gene]
})

# Subset for these
to_plot = tissue_expression_sum[Gene_ID %in% specificity_examples,]

# Scale within genes...
to_plot = to_plot[,.(tpm_scale = .SD$tpm_mean/max(.SD$tpm_mean,na.rm = T),
                     tissue),by = Gene_ID]
to_plot$TAIR_F = factor(to_plot$Gene_ID,levels = specificity_examples)
tpm_fig = ggplot(to_plot,aes(y = TAIR_F,x = tissue,fill = tpm_scale)) +
  geom_tile() +
  scale_fill_viridis() +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 10,angle = 45,hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_blank(),
        legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12,angle=45,hjust= 1)) +
  labs(fill = "TPM (Scaled)")

tau_to_plot = pergene_tau_merge[TAIR_gene %in% specificity_examples,]
tau_to_plot$gene_F = factor(tau_to_plot$TAIR_gene,levels = specificity_examples)
tau_fig = ggplot(tau_to_plot,aes(y = gene_F,x = tau_nolog)) +
  geom_bar(stat = "identity",colour = "black") +
  theme_minimal() +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.text = element_text(angle = 45,hjust = 1,size = 10)) +
  labs(x = expression(tau))

tau_demo = cowplot::plot_grid(tpm_fig,tau_fig,ncol = 2,rel_widths = c(5,1),axis = "tblr",align = "h")
tau_demo

# # Look at pleiotropy in terms of KO mutations... --------------------------
# # AT THE MOMENT THIS PARTICULAR ANALYSIS IS LIMITED BY ONLY A FEW GENES ACTUALLY HAVING KO DATA AVAILABLE...
# ko_res = list.files("data/ko")
# 
# # Run through each of these files and return the top top_N_genes associated genes...
# top_N_genes = 50
# top_ko_genes = rbindlist(lapply(ko_res,function(res){
# 
#   res_tmp = data.table(read.csv(paste0("data/ko/",res)))
# 
#   # Return top 50 genes based on pval
#   # out = data.table(Gene = res_tmp[order(Pval),SNP][1:top_N_genes],
#   #                  KO = res)
# 
#   # Return the top N genes based on var explained
#   out = data.table(Gene = res_tmp[order(-variance_explained),SNP][1:top_N_genes],
#                    KO = res)
# 
# }))
# 
# # Compile to KO counts per gene...
# KO_counts = data.table(table(top_ko_genes$Gene))[order(-N),]
# colnames(KO_counts) = c("TAIR_gene","KO_N")
# mapped_TAIR_genes = OG_map_Athal[grep("AT",OG_map_Athal$TAIR_gene),"TAIR_gene"]
# # KO_counts = rbind(KO_counts,
# #                   data.table(TAIR_gene = mapped_TAIR_genes[!mapped_TAIR_genes %in% KO_counts$TAIR_gene],
# #                              KO_N = 0))
# 
# # Merge with Orthogroup, and repeat above process to condense pvals down to OG-specific...
# KO_counts_merge = data.table(merge(OG_map_Athal[OG_map_Athal$TAIR_gene %in% mapped_TAIR_genes,c("Orthogroup","TAIR_gene")],
#                                    KO_counts,
#                                    by = "TAIR_gene"))
# KO_counts_merge$KO_pval = empPvals(KO_counts_merge$KO_N,KO_counts_merge$KO_N)
# 
# # Combine these into OG specific
# KO_OG_pvals = condenseToOGPvals(data = KO_counts_merge,
#                                 OG_col = "Orthogroup",
#                                 pval_col = "KO_pval")
# 
# # Add Z score
# KO_OG_pvals$KO_Z = qnorm(KO_OG_pvals$minP_DS_epvalue,lower.tail = F)
# 
# # Now draw the OG that we have earmarked as convergent...
# convergent_KO = KO_OG_pvals[Orthogroup %in% unique(picmin_fdr[picmin_fdr < 0.2,Orthogroup]),]
# 
# # Exclude these orthogroups, and re-draw 10000 sets
# KO_redraws = unlist(pbmclapply(1:10000,function(x){
#   set.seed(x)
#   mean(sample(KO_OG_pvals[!(Orthogroup %in% convergent_KO$Orthogroup),KO_Z],
#               nrow(convergent_KO)))
# },mc.cores = n_cores))
# 
# (mean(convergent_KO$KO_Z) - mean(KO_redraws)) / sd(KO_redraws)
# 
#### AN ALTERNATIVE WAY TO DO THIS ANALYSIS IS TO MAKE A VERSION OF TAU FOR KO USING SAME FORMULA ####
# ko_res = list.files("data/ko")
# 
# # Run through each of these files and return the top top_N_genes associated genes...
# all_ko_genes = rbindlist(pbmclapply(ko_res,function(res){
# 
#   res_tmp = data.table(read.csv(paste0("data/ko/",res)))
# 
# 
#   # Return the top N genes based on var explained
#   out = data.table(TAIR_gene = res_tmp$SNP,
#                    var_expl = res_tmp$variance_explained * 100,
#                    KO = res)
# 
# },mc.cores = n_cores))
# 
# # For all genes, calculate tau_ko
# tau_ko = all_ko_genes[, .(tau_ko = calcMankTau(.SD$var_expl),
#                           N_pheno = nrow(.SD)), by = TAIR_gene]
# 
# # Merge with OG
# tau_ko_merge = merge(tau_ko,OG_map_Athal[,c("TAIR_gene","Orthogroup")],by = "TAIR_gene")
# # Get empPvals
# tau_ko_merge$tau_ko_p = empPvals(tau_ko_merge$tau_ko,tau_ko_merge$tau_ko)
# # Condense to OG pvals
# OG_tau_ko = condenseToOGPvals(tau_ko_merge,"Orthogroup","tau_ko_p")
# OG_tau_ko$tau_Z = qnorm(OG_tau_ko$minP_DS_epvalue,lower.tail = F)
# 
# # Now draw the OG that we have earmarked as convergent...
# convergent_KO = OG_tau_ko[Orthogroup %in% unique(picmin_fdr[picmin_fdr < 0.5,Orthogroup]),]
# 
# # Exclude these orthogroups, and re-draw 10000 sets
# KO_redraws = unlist(pbmclapply(1:10000,function(x){
#   set.seed(x)
#   mean(sample(OG_tau_ko[!(Orthogroup %in% convergent_KO$Orthogroup),tau_Z],
#               nrow(convergent_KO)))
# },mc.cores = n_cores))
# 
# (mean(convergent_KO$tau_Z) - mean(KO_redraws)) / sd(KO_redraws)
# 

# Pleiotropy analysis based on GWAS hits ----------------------------------
# Sort genes on the basis of how many 'phenotypes' they are associated with...
# aragwas_snp <- fread("data/aragwas_associations_significant_permutation.csv")
# phenotypes = unique(aragwas_snp$study.phenotype.name)
# 
# # Count the number of unique associations per gene per phenotype as a proxy for pleiotropy...
# gene_unique_pheno_hits = data.table(table(unique(aragwas_snp[,.(snp.gene_name,study.phenotype.name)])[,snp.gene_name]))
# colnames(gene_unique_pheno_hits) <- c("TAIR_gene","GWAS_hits")
# gene_unique_pheno_hits = gene_unique_pheno_hits[TAIR_gene != "",]
# gene_unique_pheno_hits[order(-gene_unique_pheno_hits$GWAS_hits),]
# hist(gene_unique_pheno_hits$GWAS_hits)
# 
# # Clean gene names
# gene_unique_pheno_hits$TAIR_gene = gsub("Gene_","",gene_unique_pheno_hits$TAIR_gene)
# 
# # Count the number of hits, ignoring uniqueness...
# gene_gwas_hits = data.table(table(aragwas_snp$snp.gene_name))
# colnames(gene_gwas_hits) <- c("Gene","GWAS_hits")
# gene_gwas_hits[order(-gene_gwas_hits$GWAS_hits),]

# # If TAIR genes are missing from this list, we assume a KO_N of 0
# mapped_TAIR_genes = OG_map_Athal[grep("AT",OG_map_Athal$TAIR_gene),"TAIR_gene"]
# gene_unique_pheno_hits = rbind(gene_unique_pheno_hits,
#                                data.table(TAIR_gene = mapped_TAIR_genes[!mapped_TAIR_genes %in% gene_unique_pheno_hits$TAIR_gene],
#                                           GWAS_hits = 0))
# 
# 
# # Merge with Orthogroups and condense as previously
# gene_unique_pheno_hits_merge = merge(gene_unique_pheno_hits,OG_map_Athal[,c("TAIR_gene","Orthogroup")],by = "TAIR_gene")[order(-GWAS_hits),]
# gene_unique_pheno_hits_merge$gwas_pval = empPvals(gene_unique_pheno_hits_merge$GWAS_hits,gene_unique_pheno_hits_merge$GWAS_hits)
# OG_gwas_pvals = condenseToOGPvals(gene_unique_pheno_hits_merge,
#                                   OG_col = "Orthogroup", pval_col = "gwas_pval")
# 
# # Convert to Z-score...
# OG_gwas_pvals$gwas_Z = qnorm(OG_gwas_pvals$minP_DS_epvalue)
# 
# # Find convergent OG...
# convergent_gwas = OG_gwas_pvals[Orthogroup %in% convergent_OG,]
# 
# # Exclude these orthogroups, and re-draw 10000 sets
# GWAS_redraws = unlist(pbmclapply(1:10000,function(x){
#   set.seed(x)
#   mean(sample(OG_gwas_pvals[!(Orthogroup %in% convergent_gwas$Orthogroup),gwas_Z],
#               nrow(convergent_gwas)))
# },mc.cores = n_cores))
# 
# hist(GWAS_redraws)
# (mean(convergent_gwas$gwas_Z) - mean(GWAS_redraws)) / sd(GWAS_redraws)
# 
# 
# Pleiotropy assessed based on coexpression from ATTED-II -----------------
# Data comes from coexpression version Ath-u.c3-0
# Coexpression statistics are calculated elsewhere...

athal_node_stats = readRDS(paste0("outputs/220525_Athal_coexpression_node_stats.rds"))
centrality_vars = colnames(athal_node_stats)[-1]

# Merge with OG and remove any orthogroups that we did not test for convergence...
athal_node_stats = merge(athal_node_stats,data.table(OG_map_Athal)[,.(TAIR_gene,Orthogroup)],by = "TAIR_gene")
athal_node_stats = athal_node_stats[Orthogroup %in% unique(picmin_fdr$Orthogroup),]

# Get Z scores
centrality_Z_scores = lapply(centrality_vars,function(var){
  node_stat_vector = data.frame(athal_node_stats)[,var]
  to_test = data.table(tmp_p = empPvals(node_stat_vector,node_stat_vector),
                       Orthogroup = athal_node_stats$Orthogroup)
  
  OG_coexpressed = condenseToOGPvals(to_test,"Orthogroup","tmp_p")
  OG_coexpressed[,paste0(var,"_Z")] = qnorm(OG_coexpressed$minP_DS_epvalue,lower.tail = F)
  return(OG_coexpressed)
})
names(centrality_Z_scores) = centrality_vars

# Calculate centrality Z scores across picmin deciles
Athal_centrality_deciles = rbindlist(lapply(centrality_vars,function(var){
  tmp_merge = merge(centrality_Z_scores[[var]],picmin_deciles)
  Zcol = grep("_Z",colnames(tmp_merge))
  tmp_merge$centrality_Z = tmp_merge[,..Zcol]
  out = tmp_merge[,.(stouffersZ = stoufferZ(centrality_Z)),by = picmin_decile][order(picmin_decile)]
  out$variable = stringr::str_to_title(gsub("node_","At ",var))
  out
}))

# Merge these with our previous set...
picmin_decile_Zscores = rbind(picmin_decile_Zscores,Athal_centrality_deciles)


# Repeat climate-specific analysis for node degree ------------------------
climate_vars = unique(picmin_fdr$climate_var)
nodedegree_climate_stouffers = rbindlist(lapply(climate_vars,function(clim){
  
  convergent_OG = unique(picmin_fdr[picmin_p_adj < 0.005 & climate_var == clim,Orthogroup])
  
  # Find the tau scores associated with convergent OG
  clim_Z = centrality_Z_scores[["node_degree"]][Orthogroup %in% convergent_OG,node_degree_Z]
  clim_stouffer = sum(clim_Z)/sqrt(length(clim_Z))
  
  # Output...
  out = data.table(climate_var = clim,
                   tested_N = length(convergent_OG),
                   nodedegree_stouffer = clim_stouffer,
                   pval = pnorm(clim_stouffer,lower.tail = F))
}))
nodedegree_climate_stouffers[order(-nodedegree_stouffer),]

# Combine these with the tau climate Z and plot
nodedegree_climate_stouffers$stouffer = nodedegree_climate_stouffers$nodedegree_stouffer
tau_climate_stouffers$stouffer = tau_climate_stouffers$tau_stouffer
plot_climate_Z = rbind(nodedegree_climate_stouffers[,.(climate_var,stouffer)],
                       tau_climate_stouffers[,.(climate_var,stouffer)])
plot_climate_Z$variable = rep(c("At Node Degree","Tau (Specificity)"),each = nrow(plot_climate_Z)/2)
plot_climate_Z$climate_F = factor(stringr::str_to_title(gsub("_"," ",plot_climate_Z$climate_var)),
                                  levels = stringr::str_to_title(gsub("_"," ",tau_climate_stouffers[order(stouffer),climate_var])))

climate_specific_Z = ggplot(plot_climate_Z,aes(y = climate_F,x = stouffer,fill = variable)) +
  geom_bar(stat = 'identity',position = 'dodge') +
  geom_vline(xintercept = qnorm(0.025,lower.tail = F),linetype = 'dashed') +
  scale_fill_manual(values = c("Tau (Specificity)" = "#E69A8DFF",
                               "At Node Degree" = "#5F4B8BFF")) +
  theme_minimal() +
  theme(panel.grid.minor.y = element_blank(),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = "top",
        legend.title = element_blank(),
        legend.direction = "vertical") +
  labs(x = "Stouffer's Z\n(Repeatable Orthogroups)")
climate_specific_Z
# # Fetch mean Z score for each centrality measure
# convergent_centrality_Z = sapply(centrality_vars,function(var) {
#   mean(unlist(centrality_Z_scores[[var]][Orthogroup %in% convergent_OG,paste0(var,"_Z"),with=F]))
# })
# convergent_centrality_Z
# 
# # Exclude these orthogroups, and re-draw 10000 sets
# coexpress_redraws = rbindlist(pbmclapply(1:10000,function(x){
#   set.seed(x)
#   # mean(sample(OG_coexpressed[!(Orthogroup %in% convergent_OG),coexpressed_Z],
#   #             nrow(convergent_coexpressed)))
#   
#   random_centrality_Z = sapply(centrality_vars,function(var) {
#     mean(sample(unlist(centrality_Z_scores[[var]][!(Orthogroup %in% convergent_OG),paste0(var,"_Z"),with=F]),
#                 size = length(convergent_OG)))
#   })
#   
#   out = data.table(random_centrality_Z)
#   out$centrality = centrality_vars
#   out
# },mc.cores = n_cores))
# 
# # Calculate final Z scores and plot results...
# convergent_coexpress_Zs = sapply(centrality_vars,function(var){
#   (convergent_centrality_Z[var] - mean(coexpress_redraws[centrality == var,random_centrality_Z])) / sd(coexpress_redraws[centrality == var,random_centrality_Z])
# })
# convergent_coexpress_Zs
# 
# # Get pvals
# convergent_coexpress_pvals = data.table(centrality = paste0(stringr::str_to_title(gsub("_"," ",centrality_vars))," - Arabidopsis"),
#                                         epval = pnorm(convergent_coexpress_Zs,lower.tail = F),
#                                         Z = convergent_centrality_Z)
# 
# coexpress_redraws$centrality = paste0(stringr::str_to_title(gsub("_"," ",coexpress_redraws$centrality))," - Arabidopsis")
# 
# ggplot(data = coexpress_redraws,aes(random_centrality_Z))+
#   geom_histogram(bins = 20) +
#   # geom_vline(xintercept = convergent_coexpress_pvals$epval,colour = "red2",size = 2) +
#   geom_vline(data = convergent_coexpress_pvals,aes(xintercept = Z),colour = "red2",size = 2) +
#   facet_wrap(~centrality) + 
#   theme_minimal() +
#   theme(axis.title = element_text(size = 18),
#         axis.text = element_text(size = 14),
#         title = element_text(size = 20),
#         strip.text = element_text(size = 24)) + 
#   labs(x = expression(Mean~Coexpression~Z),y = "Count")
# 
# # Add these re-draws onto our redraw data
# colnames(coexpress_redraws) = c("randomZ","variable")
# colnames(convergent_coexpress_pvals) = c("variable","epval","obsZ")
# re_draws_to_plot = rbind(re_draws_to_plot,coexpress_redraws)
# convergent_means = rbind(convergent_means,convergent_coexpress_pvals[,.(obsZ,variable)])
# 
# Merge all the Z scores together for exporting
OG_arabidopsis_coexpress_Z = data.table(Orthogroup = centrality_Z_scores$node_betweenness$Orthogroup,
                                        At_node_betweenness_Z = centrality_Z_scores$node_betweenness$node_betweenness_Z,
                                        At_node_strength_Z = centrality_Z_scores$node_strength$node_strength_Z,
                                        At_node_degree_Z = centrality_Z_scores$node_degree$node_degree_Z,
                                        At_node_closeness_Z = centrality_Z_scores$node_closeness$node_closeness_Z)


# What is effect size difference of node degree ---------------------------
convergent_OG = unique(picmin_fdr[picmin_p_adj < 0.005,Orthogroup])

OG_max_degree = athal_node_stats[,.(max_degree = max(node_degree),
                                    mean_degree = mean(node_degree),
                                    max_strength = max(node_strength),
                                    mean_strength = mean(node_strength)),by = Orthogroup]

# How do these differ?
mean(OG_max_degree[Orthogroup %in% convergent_OG,max_degree]) / mean(OG_max_degree[!Orthogroup %in% convergent_OG,max_degree])
mean(OG_max_degree[Orthogroup %in% convergent_OG,mean_degree]) / mean(OG_max_degree[!Orthogroup %in% convergent_OG,mean_degree])

mean(OG_max_degree[Orthogroup %in% convergent_OG,max_strength]) / mean(OG_max_degree[!Orthogroup %in% convergent_OG,max_strength])
mean(OG_max_degree[Orthogroup %in% convergent_OG,mean_strength]) / mean(OG_max_degree[!Orthogroup %in% convergent_OG,mean_strength])

# And plot deciles
OG_degree = OG_arabidopsis_coexpress_Z
OG_degree$convergent = "No"
OG_degree[Orthogroup %in% convergent_OG,convergent := "Yes"]

# What are the deciles...
OG_degree$deciles = cut_number(OG_degree$At_node_degree_Z,10)

# Plot these...
degree_deciles = data.table(obs = table(OG_degree$deciles) / nrow(OG_degree),
                            convergent = table(OG_degree[convergent == "Yes",deciles]) / nrow(OG_degree[convergent == "Yes",]))
degree_deciles$decile = factor(1:10,levels = 1:10)

degree_decile_plot = melt(degree_deciles[,.(decile,obs.N,convergent.N)]) |> 
  ggplot(aes(y = value, x = factor(decile),fill = variable)) + 
  geom_bar(stat = 'identity',position = 'dodge') +
  theme_minimal() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        title = element_text(size = 16)) +
  labs(y = "Proportion",x = "Decile",fill = "Repeated?") +
  scale_fill_discrete(breaks = c("obs.N","convergent.N"),
                      labels = c("Expected","Repeated")) +
  ggtitle("Arabidopsis Node Degree")
degree_decile_plot
# Plot each as subfig
pdf("figs/FigureSX_tau_degree_deciles.pdf",width = 10,height = 8)
cowplot::plot_grid(tau_decile_plot,
                   degree_decile_plot,
                   ncol = 1,labels = "AUTO",label_size = 32,align = 'v',axis = 'tblr')
dev.off()



# Pleiotropy ATTED-II using Medicago --------------------------------------
# Data comes from coexpression version Mtr-u.c3-0
# Coexpression statistics are calculated elsewhere...

mtrunc_node_stats = readRDS(paste0("outputs/220927_Mtrunc_coexpression_node_stats.rds"))
centrality_vars = colnames(mtrunc_node_stats)[-1]

# Merge with OG and remove any orthogroups that we did not test for convergence...
mtrunc_node_stats = merge(mtrunc_node_stats,data.table(OG_map_Mtrunc)[,.(biomart_gene,Orthogroup)],by = "biomart_gene")
mtrunc_node_stats = mtrunc_node_stats[Orthogroup %in% unique(picmin_fdr$Orthogroup),]

# Get Z scores
centrality_Z_scores = lapply(centrality_vars,function(var){
  node_stat_vector = data.frame(mtrunc_node_stats)[,var]
  to_test = data.table(tmp_p = empPvals(node_stat_vector,node_stat_vector),
                       Orthogroup = mtrunc_node_stats$Orthogroup)
  
  OG_coexpressed = condenseToOGPvals(to_test,"Orthogroup","tmp_p")
  OG_coexpressed[,paste0(var,"_Z")] = qnorm(OG_coexpressed$minP_DS_epvalue,lower.tail = F)
  return(OG_coexpressed)
})
names(centrality_Z_scores) = centrality_vars

# Calculate centrality Z scores across picmin deciles
Mtrunc_centrality_deciles = rbindlist(lapply(centrality_vars,function(var){
  tmp_merge = merge(centrality_Z_scores[[var]],picmin_deciles)
  Zcol = grep("_Z",colnames(tmp_merge))
  tmp_merge$centrality_Z = tmp_merge[,..Zcol]
  out = tmp_merge[,.(stouffersZ = stoufferZ(centrality_Z)),by = picmin_decile][order(picmin_decile)]
  out$variable = stringr::str_to_title(gsub("node_","Mt ",var))
  out
}))

# Merge these with our previous set...
picmin_decile_Zscores = rbind(picmin_decile_Zscores,Mtrunc_centrality_deciles)


# # Fetch mean Z score for each centrality measure
# convergent_centrality_Z = sapply(centrality_vars,function(var) {
#   mean(unlist(centrality_Z_scores[[var]][Orthogroup %in% convergent_OG,paste0(var,"_Z"),with=F]))
# })
# convergent_centrality_Z
# 
# # Exclude these orthogroups, and re-draw 10000 sets
# coexpress_redraws = rbindlist(pbmclapply(1:10000,function(x){
#   set.seed(x)
#   # mean(sample(OG_coexpressed[!(Orthogroup %in% convergent_OG),coexpressed_Z],
#   #             nrow(convergent_coexpressed)))
#   
#   random_centrality_Z = sapply(centrality_vars,function(var) {
#     mean(sample(unlist(centrality_Z_scores[[var]][!(Orthogroup %in% convergent_OG),paste0(var,"_Z"),with=F]),
#                 size = length(convergent_OG)))
#   })
#   
#   out = data.table(random_centrality_Z)
#   out$centrality = centrality_vars
#   out
# },mc.cores = n_cores))
# 
# # Calculate final Z scores and plot results...
# convergent_coexpress_Zs = sapply(centrality_vars,function(var){
#   (convergent_centrality_Z[var] - mean(coexpress_redraws[centrality == var,random_centrality_Z])) / sd(coexpress_redraws[centrality == var,random_centrality_Z])
# })
# convergent_coexpress_Zs
# 
# # Get pvals
# convergent_coexpress_pvals = data.table(centrality = paste0(stringr::str_to_title(gsub("_"," ",centrality_vars))," - Medicago"),
#                                         epval = pnorm(convergent_coexpress_Zs,lower.tail = F),
#                                         Z = convergent_centrality_Z)
# convergent_coexpress_pvals
# coexpress_redraws$centrality = paste0(stringr::str_to_title(gsub("_"," ",coexpress_redraws$centrality))," - Medicago")
# 
# ggplot(data = coexpress_redraws,aes(random_centrality_Z))+
#   geom_histogram(bins = 20) +
#   # geom_vline(xintercept = convergent_coexpress_pvals$epval,colour = "red2",size = 2) +
#   geom_vline(data = convergent_coexpress_pvals,aes(xintercept = Z),colour = "red2",size = 2) +
#   facet_wrap(~centrality) + 
#   theme_minimal() +
#   theme(axis.title = element_text(size = 18),
#         axis.text = element_text(size = 14),
#         title = element_text(size = 20),
#         strip.text = element_text(size = 24)) + 
#   labs(x = expression(Mean~Coexpression~Z),y = "Count")
# 
# # Add these re-draws onto our redraw data
# colnames(coexpress_redraws) = c("randomZ","variable")
# colnames(convergent_coexpress_pvals) = c("variable","epval","obsZ")
# re_draws_to_plot = rbind(re_draws_to_plot,coexpress_redraws)
# convergent_means = rbind(convergent_means,convergent_coexpress_pvals[,.(obsZ,variable)])
# convergent_means
# 
# Merge all the Z scores together for exporting
OG_medicago_coexpress_Z = data.table(Orthogroup = centrality_Z_scores$node_betweenness$Orthogroup,
                                     Mt_node_betweenness_Z = centrality_Z_scores$node_betweenness$node_betweenness_Z,
                                     Mt_node_strength_Z = centrality_Z_scores$node_strength$node_strength_Z,
                                     Mt_node_degree_Z = centrality_Z_scores$node_degree$node_degree_Z,
                                     Mt_node_closeness_Z = centrality_Z_scores$node_closeness$node_closeness_Z)

# Combine all of the Z scores by OG...
OG_pleiotropy_Z = merge(OG_tau1[,.(Orthogroup,tau_Z)],OG_medicago_coexpress_Z,by = "Orthogroup")
OG_pleiotropy_Z = merge(OG_pleiotropy_Z,OG_arabidopsis_coexpress_Z,by = "Orthogroup")
OG_pleiotropy_Z

# # Plot pleiotropy results... ----------------------------------------------
# re_draws_to_plot$variable = gsub("tau","Tissue Specificity - Tau",re_draws_to_plot$variable)
# convergent_means$variable = gsub("tau","Tissue Specificity - Tau",convergent_means$variable)
# 
# re_draws_to_plot$var_F = factor(re_draws_to_plot$variable,levels = unique(re_draws_to_plot$variable))
# convergent_means$var_F = factor(convergent_means$variable,levels = unique(convergent_means$variable))
# 
# # Calculate means and upper/lower quantiles...
# re_draws_means_quants = re_draws_to_plot[,.(meanZ = mean(randomZ),
#                                             q5 = quantile(randomZ,probs = 0.05),
#                                             q95 = quantile(randomZ,probs = 0.95)),by = var_F]
#
# ggplot(re_draws_to_plot,aes(y = var_F, x = randomZ)) +
#   geom_density_ridges() +
#   geom_point(data = convergent_means,aes(y = var_F, x = obsZ),colour = "forestgreen",size = 5) +
#   scale_y_discrete(limits=rev)
# 
# # Add some fairly arbitrary pval
# convergent_Z$pval = pnorm(convergent_Z$stoufferZ,lower.tail = F) * 2
# convergent_Z$var_F = factor(convergent_Z$variable,levels = convergent_Z$variable)
# 
# # Add in asterices
# convergent_Z$pval_tag = ""
# convergent_Z[pval < 0.05,pval_tag := "*"]
# convergent_Z[pval < 0.01,pval_tag := "**"]
# convergent_Z[pval < 0.001,pval_tag := "***"]
# 
# pleiotropy_results_fig = ggplot(convergent_Z,aes(y = var_F,x = stoufferZ)) + 
#   geom_bar(stat = 'identity',colour = 'black',fill = 'gold2') +
#   scale_y_discrete(limits=rev) +
#   theme_minimal() + 
#   theme(axis.title.y = element_blank(),
#         axis.title.x = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         panel.grid.major.y = element_blank()) +
#   labs(x = "Stouffer's Z\n(Repeatable Orthogroups)") +
#   geom_text(data = convergent_Z,aes(x = convergent_Z$stoufferZ + 1,label = pval_tag,y = var_F),size = 10) +
#   # xlim(-0.15,0.35) +
#   geom_vline(xintercept = qnorm(0.025,lower.tail = F),linetype = 'dashed')
# pleiotropy_results_fig

# pleiotropy_results_fig = ggplot(re_draws_means_quants,aes(y = var_F,x = meanZ)) + 
#   geom_point(data = convergent_means,aes(y = var_F, x = obsZ),colour = "forestgreen",size = 5,shape = 17) +
#   geom_point(size = 5) +
#   geom_segment(aes(x = q5,xend = q95,y = var_F,yend = var_F),size = 1.5) +
#   scale_y_discrete(limits=rev) +
#   theme_minimal() + 
#   theme(axis.title.y = element_blank(),
#         axis.title.x = element_text(size = 14),
#         axis.text = element_text(size = 12),
#         panel.grid.major.y = element_blank()) +
#   labs(x = "Mean Orthogroup Z score") +
#   geom_text(data = variable_pvals,aes(x = 0.3,label = pval_tag,y = var_F),size = 10) +
#   xlim(-0.15,0.35) +
#   geom_vline(xintercept = 0,linetype = "dotted")

# Pleiotropy results figure based on the picmin deciles
picmin_decile_Zscores$variable_labs = gsub(" ","\n",picmin_decile_Zscores$variable)
picmin_decile_Zscores$Z_col <- "None"
picmin_decile_Zscores[stouffersZ > qnorm(1 - (0.05 / 2)),Z_col := "Pos"]
picmin_decile_Zscores[stouffersZ < -qnorm(1 - (0.05 / 2)),Z_col := "Neg"]

pleiotropy_results_fig = ggplot(picmin_decile_Zscores,aes(y = picmin_decile,x = stouffersZ,
                                 size = abs(stouffersZ),colour = Z_col)) +
  geom_point() +
  facet_wrap(~variable_labs,nrow = 1) +
  scale_y_discrete(limits = rev) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = c(qnorm(1 - (0.05 / 2)),-qnorm(1 - (0.05 / 2))),
             linetype = 'dotted') +
  scale_colour_manual(values = c("None" = "black",
                                 "Pos" = "forestgreen",
                                 "Neg" = "red4")) +
  labs(y = "Evidence for repeatability\n(min PicMin p-value)",x = "Stouffer's Z") +
  theme_bw() +
  theme(legend.position = 'none',
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        strip.text = element_text(size = 12),
        strip.background = element_blank())

# Combine with other figs...
source("R/plot_network_centrality_demo.R")

# pdf("figs/Figure4_pleiotropy_and_convergence_results.pdf",width = 14,height = 10)
# cowplot::plot_grid(
#   cowplot::plot_grid(tau_demo,
#                      centrality_combined,
#                      ncol=1, labels = c("A","B"),label_size = 32,label_y = c(1,1.2)),
#   cowplot::plot_grid(pleiotropy_results_fig,NULL,nrow = 2),
#   ncol = 2,
#   labels = c("","C"),label_size = 32
# )
# dev.off()

# Big export for plotting elsewhere...
saveRDS(list(OG_pleiotropy_Z = OG_pleiotropy_Z,
             tau_demo = tau_demo,
             centrality_combined = centrality_combined,
             pleiotropy_results_fig = pleiotropy_results_fig,
             climate_specific_Z = climate_specific_Z),
        paste0('outputs/',output_name,'_all_pleiotropy_figs_and_OG_Zscores.rds'))


# # Repeat analysis for adaptation OG
# nonconvergent_coexpressed = OG_coexpressed[Orthogroup %in% nonconvergent_OG,]
# hist(nonconvergent_coexpressed$coexpressed_Z)
# 
# # Exclude these orthogroups, and re-draw 10000 sets
# coexpress_redraws2 = unlist(pbmclapply(1:10000,function(x){
#   set.seed(x)
#   mean(sample(OG_coexpressed[!(Orthogroup %in% nonconvergent_OG),coexpressed_Z],
#               nrow(nonconvergent_coexpressed)))
# },mc.cores = n_cores))
# 
# # Calculate final Z score and plot results...
# hist(coexpress_redraws2)
# nonconvergent_coexpress_Z = (mean(nonconvergent_coexpressed$coexpressed_Z) - mean(coexpress_redraws2)) / sd(coexpress_redraws2)
# nonconvergent_coexpress_Z
# 
# data.table(coexpress_redraws) |>
#   ggplot(aes(coexpress_redraws))+
#   geom_histogram(bins = 20) +
#   geom_vline(xintercept = mean(convergent_coexpressed$coexpressed_Z),colour = "red2",size = 2) +
#   geom_vline(xintercept = mean(nonconvergent_coexpressed$coexpressed_Z),colour = "gold2",size = 2) +
#   theme_minimal() +
#   theme(axis.title = element_text(size = 18),
#         axis.text = element_text(size = 14),
#         title = element_text(size = 20)) + 
#   labs(x = expression(Mean~Coexpression~Z),y = "Count") +
#   ggtitle(paste0("Convergent OG p = ",round(pnorm(convergent_coexpress_Z,lower.tail = F),4)))

# # Analysis of gene duplication and paralogy on convergence ----------------
# 
# # # We can use the OG map to infer the rate of paralogy
# # OG_map = OG_map[genome != "Pmenziesii",]
# # OG_paralog_N = OG_map[,.(paralogN = nrow(.SD)), by = .(Orthogroup,genome)]
# # OG_paralog_stats = OG_paralog_N[paralogN <= 10,.(max_paralogN = max(.SD$paralogN),
# #                                                  mean_paralogN = mean(.SD$paralogN),
# #                                                  total_OG_size = sum(.SD$paralogN),
# #                                                  genomeN = nrow(.SD),
# #                                                  # median_paralogN = median(.SD$paralogN),
# #                                                  duplicatedN = sum(.SD$paralogN > 1)),by=Orthogroup]
# # Infer paralog rates on the basis of the paralogs we tested for convergence...
# OG_pvals_tested = OG_pvals[Orthogroup %in% unique(picmin_fdr$Orthogroup),]
# OG_paralog_stats = OG_pvals_tested[climate == "mean_temp",][,.(max_paralogN = max(.SD$Ngenes_per_species/.SD$Ndatasets_per_species),
#                                                               mean_paralogN = mean(.SD$Ngenes_per_species/.SD$Ndatasets_per_species),
#                                                               total_OG_size = sum(.SD$Ngenes_per_species/.SD$Ndatasets_per_species),
#                                                               genomeN = nrow(.SD),
#                                                               # median_paralogN = median(.SD$paralogN),
#                                                               duplicatedN = sum((.SD$Ngenes_per_species/.SD$Ndatasets_per_species) > 1)),
#                                                             by=Orthogroup]
# 
# # Filter this for orthogroups that we actually explored for convergence...
# OG_paralog_stats = OG_paralog_stats[Orthogroup %in% unique(picmin_fdr$Orthogroup),]
# 
# # # Merge these stats with the picmin_fdr and plot...
# picmin_fdr_paralogs = merge(picmin_fdr,OG_paralog_stats,by = "Orthogroup")
# # ggplot(picmin_fdr_paralogs,aes(max_paralogN,-log10(p)))+
# #   geom_point()+
# #   stat_density2d_filled(alpha=0.5)
# # ggplot(picmin_fdr_paralogs,aes((mean_paralogN),(-log10(p))))+
# #   geom_point()+
# #   stat_density2d_filled(alpha=0.5)
# # ggplot(picmin_fdr_paralogs,aes(duplicatedN,(-log10(p))))+
# #   geom_point()+
# #   stat_density2d_filled(alpha=0.5)
# 
# # Fetch stats for our convergent_OG
# convergent_OG = unique(picmin_fdr[picmin_fdr < 0.2,Orthogroup])
# # convergent_OG = unique(convergent_OG[duplicated(convergent_OG)])
# convergent_OG_paralog_stats = unique(picmin_fdr_paralogs[Orthogroup %in% convergent_OG,.(Orthogroup,max_paralogN,mean_paralogN,duplicatedN)])
# 
# # Exclude these orthogroups, and re-draw 10000 sets
# to_draw = unique(picmin_fdr_paralogs$Orthogroup[!picmin_fdr_paralogs$Orthogroup %in% convergent_OG])
# paralog_redraws = rbindlist(pbmclapply(1:10000,function(x){
#   set.seed(x)
#   
#   tmp = unique(picmin_fdr_paralogs[Orthogroup %in% sample(to_draw,length(convergent_OG)),
#                                    .(Orthogroup,max_paralogN,mean_paralogN,duplicatedN)])
#   
#   data.table(`Max Paralog N` = mean(tmp$max_paralogN),
#              `Mean Paralog N` = mean(tmp$mean_paralogN),
#              `Duplicated N` = mean(tmp$duplicatedN),
#              draw = x)
# },mc.cores = n_cores))
# 
# # Test each one...
# # Max Paralog N
# empPvals(mean(convergent_OG_paralog_stats$max_paralogN),stat0 = paralog_redraws$mean_maxparalogN)
# # Mean Paralog N
# empPvals(mean(convergent_OG_paralog_stats$mean_paralogN),stat0 = paralog_redraws$mean_meanparalogN)
# # Number of duplications Duplications
# empPvals(mean(convergent_OG_paralog_stats$duplicatedN),stat0 = paralog_redraws$mean_duplicateN)
# 
# # Plot verticals
# paralog_stats_to_plot = data.table(mean_val = c(mean(convergent_OG_paralog_stats$max_paralogN),
#                                                 mean(convergent_OG_paralog_stats$mean_paralogN),
#                                                 mean(convergent_OG_paralog_stats$duplicatedN)),
#                                    variable = c("Max Paralog N",
#                                                 "Mean Paralog N",
#                                                 "Duplicated N"))
# 
# # Plot duplication results...
# melt(paralog_redraws)[variable != "draw",] |>
#   ggplot(aes(value)) +
#   geom_histogram() +
#   facet_wrap(~variable,scales = "free") +
#   geom_vline(data = paralog_stats_to_plot,aes(xintercept = mean_val),colour="red2",size = 2) +
#   theme_minimal() +
#   theme(axis.title = element_text(size = 18),
#         axis.text = element_text(size = 14),
#         title = element_text(size = 20),
#         strip.text = element_text(size = 24)) + 
#   labs(x = "Mean Paralog Stats",y = "Count")
# 
# # We also want to see whether the probability of a gene duplication is associated with specific datasets contributing to convergence....
# OG_pvals_tested$paralogN = OG_pvals_tested$Ngenes_per_species/OG_pvals_tested$Ndatasets_per_species
# OG_pvals_tested$adaptive = "No"
# OG_pvals_tested[max_sdP_DS < 0.1, adaptive := "Yes"]
# 
# OG_pvals_tested[,median(paralogN),by = adaptive]
# 
# # First keep only convergent sets, and assess diff in paralogs between adaptive/not adaptive within convergent
# OG_pvals_tested$OG_climate = paste0(OG_pvals_tested$Orthogroup,"-",OG_pvals_tested$climate)
# to_test = apply(picmin_fdr[picmin_fdr < 0.2,.(Orthogroup,climate_var)],1,paste,collapse= "-")
# OG_pvals_tested$is_duplicated = "No"
# OG_pvals_tested[paralogN > 1, is_duplicated := "Yes"]
# 
# OG_pvals_convergent = OG_pvals_tested[OG_climate %in% to_test,]
# 
# # Average paralogN
# OG_pvals_convergent[,mean(paralogN),by = adaptive]
# 
# # Likelihood of duplication...
# plot(table(OG_pvals_convergent[,.(adaptive,is_duplicated)]))
# plot(table(OG_pvals_tested[,.(adaptive,is_duplicated)]))
# 
# # Possibly redundant section with LMM of duplication ----------------------


# # First assign each species in picmin_fdr with its relevant genome...
# OG_paralog_N$genome = gsub("Atubercatus","Atuberculatus",OG_paralog_N$genome)
# OG_pvals_genome = merge(OG_pvals,species_genome_map,by = "species")
# OG_pvals_genome = merge(OG_pvals_genome,OG_paralog_N,by = c("genome","Orthogroup"))
# 
# # Now extract convergent orthogroups and ask whether those species contributing to single show greater duplication...
# fdr_thresh = 0.5
# OG_pvals_genome$OG_climate = paste0(OG_pvals_genome$Orthogroup,"-",OG_pvals_genome$climate)
# picmin_fdr$OG_climate = paste0(picmin_fdr$Orthogroup,"-",picmin_fdr$climate_var)
# OG_pvals_genome_convergent = OG_pvals_genome[OG_climate %in% picmin_fdr[picmin_fdr < fdr_thresh,OG_climate],]
# 
# # Group species into either contributing, or not contributing...
# OG_pvals_genome_convergent$convergent = "Not-convergent"
# OG_pvals_genome_convergent[max_sdP_DS < 0.1,convergent := "Convergent"]
# ggplot(OG_pvals_genome_convergent,aes(x = paralogN, fill = convergent)) +
#   geom_density()
# ggplot(OG_pvals_genome_convergent,aes(y = log(paralogN), x = convergent)) +
#   geom_boxplot()
# 
# # Model this according to a poisson glmm
# # Fit model
# paralog_mm = lme4:::glmer(paralogN ~ convergent + (1 | species),data = OG_pvals_genome_convergent,family = poisson(link = "log"))
# 
# # Check model
# paralog_model_check = performance::check_model(paralog_mm)
# paralog_model_check
# library(DHARMa)
# testDispersion(paralog_mm)
# simulationOutput <- simulateResiduals(fittedModel = paralog_mm, plot = F)
# residuals(simulationOutput)
# residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
# plot(simulationOutput)
# 
# # Signif overdispersion, switch to NB model
# OG_pvals_genome_convergent$`Paralog N` = OG_pvals_genome_convergent$paralogN
# 
# paralog_mm2 = lme4:::glmer.nb(paralogN ~ convergent + (1 | species),data = OG_pvals_genome_convergent)
# summary(paralog_mm2)
# 
# tab_model(paralog_mm2)
# plot(effects::allEffects(paralog_mm2))
# 
# # Also simpler binomial model of duplicated vs not-duplicated...
# OG_pvals_genome_convergent$duplicated = "No"
# OG_pvals_genome_convergent[paralogN > 1, duplicated := "Yes"]
# paralog_bin_mod = lme4:::glmer(factor(duplicated) ~ convergent + (1 | species),data = OG_pvals_genome_convergent,family = binomial)
# performance::check_model(paralog_bin_mod)
# plot(effects::allEffects(paralog_bin_mod))
# 
# # How does this compare across all OG...
# OG_pvals_genome$adaptive = "No"
# OG_pvals_genome[max_sdP_DS < 0.1, adaptive := "Yes"]
# OG_pvals_genome$duplicated = "No"
# OG_pvals_genome[paralogN > 1, duplicated := "Yes"]
# 
# paralog_bin_mod2 = lme4:::glmer(factor(duplicated) ~ adaptive + (1 | species),data = OG_pvals_genome,family = binomial)
# performance::check_model(paralog_bin_mod2)
# plot(effects::allEffects(paralog_bin_mod2))
# 
# # performance::check_model(paralog_mm2)
# # testDispersion(paralog_mm2)
# # simulationOutput <- simulateResiduals(fittedModel = paralog_mm2, plot = F)
# # residuals(simulationOutput)
# # residuals(simulationOutput, quantileFunction = qnorm, outlierValues = c(-7,7))
# # plot(simulationOutput)

# # Orthogroup Functions through KEGG ---------------------------------------
# library(clusterProfiler)
# 
# # Fetch the KEGG database...
# Ath_KEGG = download_KEGG('ath')
# 
# # Make some fixes to OG_map_Athal and merge to the Ath_KEGG map...
# OG_map_Athal[grep("Arth",OG_map_Athal$gene_name_noGenome),"TAIR_gene"] <- gsub('gene-','',OG_map_Athal[grep("Arth",OG_map_Athal$gene_name_noGenome),"gene_name_noGenome"])
# kegg_genes = Ath_KEGG$KEGGPATHID2EXTID
# colnames(kegg_genes) = c("kegg_name","TAIR_gene")
# kegg_genes = data.table(merge(kegg_genes,OG_map_Athal[,c("Orthogroup","TAIR_gene")]))
# 
# # Set up kegg_paths as a list object which describes the orthogroups associated with each 
# kegg_paths = data.table(Ath_KEGG$KEGGPATHID2NAME)
# 
# # Keep only unique Orthogroup-KEGG pairs...
# kegg_OG = unique(kegg_genes[,.(Orthogroup,kegg_name)])
# kegg_list = lapply(unique(kegg_OG$kegg_name),function(path) kegg_OG[kegg_name == path,Orthogroup])
# names(kegg_list) = sapply(unique(kegg_OG$kegg_name),function(path) kegg_paths[from == path,to])
# 
# # What KEGG paths are associated with convergent OG?
# fdr_cutoff = 0.5
# 
# # Run this for all individual climate variables analysed...
# convergent_OG_list = lapply(unique(picmin_fdr$climate_var),function(var) picmin_fdr[picmin_fdr < fdr_cutoff &
#                                                                                       climate_var == var,Orthogroup])
# 
# kegg_res_list = lapply(convergent_OG_list,function(convergent_OG){
#   # print(length(convergent_OG))
#   if(length(convergent_OG) > 0){
#     
#     # Limit to KEGG paths we can actually look at...
#     kegg_to_test = kegg_paths[from %in% kegg_OG[Orthogroup %in% convergent_OG,kegg_name],to]
#     
#     if(length(kegg_to_test) > 0){
#       convergent_OG_enrich = bc3net::enrichment(genes = convergent_OG,
#                                                 reference = unique(picmin_fdr$Orthogroup),
#                                                 genesets = kegg_list[names(kegg_list) %in% kegg_paths[from %in% kegg_OG[Orthogroup %in% convergent_OG,kegg_name],to]],
#                                                 verbose = F)
#     }
#   } 
#   # head(convergent_OG_enrich)
# })
# names(kegg_res_list) = unique(picmin_fdr$climate_var)
# 

# #### SOMEWHERE ELSE FROM WHICH WE CAN PULL KEGG ANNOTATIONS IS OBVS BIOMART
# ensembl <- useMart(biomart = "plants_mart",host="https://plants.ensembl.org")
# athal_biomart <- useDataset("athaliana_eg_gene",mart=ensembl)
# athal_universe <- getBM(attributes = c("tair_locus","ensembl_gene_id","chromosome_name","start_position","end_position",
#                                        "hmmpanther","go_id","plant_reactome_pathway","plant_reactome_reaction"),
#                         mart=athal_biomart)
# 
# # Merge with Orthogroups...
# tmp = data.table(OG_map_Athal)
# tmp$ensembl_gene_id = tmp$TAIR_gene
# athal_universe_OG = data.table(merge(athal_universe,tmp[,.(ensembl_gene_id,Orthogroup)]))
# 
# #### Panther enrichment...
# panther_OG = unique(athal_universe_OG[,.(hmmpanther,Orthogroup)])
# # panther_OG = separate(panther_OG,hmmpanther,sep = ':',into = c("hmmpanther","extra"))
# # panther_OG = unique(panther_OG[,.(hmmpanther,Orthogroup)])
# 
# # Make list
# panther_list = lapply(unique(panther_OG$hmmpanther),function(x) panther_OG[hmmpanther == x,Orthogroup])
# names(panther_list) = unique(panther_OG$hmmpanther)
# 
# # What PANTHER paths are associated with convergent OG?
# fdr_cutoff = 0.5
# convergent_OG_list = lapply(unique(picmin_fdr$climate_var),function(var) picmin_fdr[picmin_fdr < fdr_cutoff &
#                                                                                       climate_var == var,Orthogroup])
# 
# panther_res_list = lapply(convergent_OG_list,function(convergent_OG){
#   # print(length(convergent_OG))
#   if(length(convergent_OG) > 0){
#     
#     # Limit to KEGG paths we can actually look at...
#     panther_to_test = unique(panther_OG[Orthogroup %in% convergent_OG,hmmpanther])
#     
#     if(length(panther_to_test) > 0){
#       convergent_OG_enrich = bc3net::enrichment(genes = convergent_OG,
#                                                 reference = unique(picmin_fdr$Orthogroup),
#                                                 genesets = panther_list[names(panther_list) %in% panther_to_test],
#                                                 verbose = F)
#       
#     }
#   }
# })
# 
# 
# #### GO enrichment...
# GO_OG = unique(athal_universe_OG[,.(go_id,Orthogroup)])
# 
# # Make list
# GO_list = pbmclapply(unique(GO_OG$go_id),function(x) GO_OG[go_id == x,Orthogroup],mc.cores = 4)
# names(GO_list) = unique(GO_OG$go_id)
# 
# # What GO paths are associated with convergent OG?
# fdr_cutoff = 0.5
# convergent_OG_list = lapply(unique(picmin_fdr$climate_var),function(var) picmin_fdr[picmin_fdr < fdr_cutoff &
#                                                                                       climate_var == var,Orthogroup])
# 
# GO_res_list = lapply(convergent_OG_list,function(convergent_OG){
#   # print(length(convergent_OG))
#   if(length(convergent_OG) > 0){
#     
#     # Limit to KEGG paths we can actually look at...
#     GO_to_test = unique(GO_OG[Orthogroup %in% convergent_OG,go_id])
#     GO_to_test = GO_to_test[GO_to_test != ""]
#     
#     if(length(GO_to_test) > 0){
#       convergent_OG_enrich = bc3net::enrichment(genes = convergent_OG,
#                                                 reference = unique(picmin_fdr$Orthogroup),
#                                                 genesets = GO_list[names(GO_list) %in% GO_to_test],
#                                                 verbose = F)
#       #Only return fdr < 0.1
#       convergent_OG_enrich[convergent_OG_enrich$padj < 0.1,]
#     }
#   }
# })



