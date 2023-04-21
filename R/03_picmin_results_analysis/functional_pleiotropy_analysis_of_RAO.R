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

# Merge with the Orthogroups from OG map
pergene_tau_merge = merge(data.table(OG_map_Athal)[,.(TAIR_gene,Orthogroup)],tissue_expression_tau,by="TAIR_gene")
# Filter for Orthogroups that we actually test for convergence...
pergene_tau_merge = pergene_tau_merge[Orthogroup %in% unique(picmin_fdr$Orthogroup),]

# Convert tau scores into empirical pvals then Z scores, correcting through DS for Orthogroup...
# Looking for low tau...
pergene_tau_merge$tau_p1 = empPvals(-pergene_tau_merge$tau_nolog,-pergene_tau_merge$tau_nolog)
OG_tau1 = condenseToOGPvals(pergene_tau_merge,"Orthogroup","tau_p1")
OG_tau1$tau_Z = qnorm(OG_tau1$minP_DS_epvalue,lower.tail = F)

# Take most significant PicMin Orthogroups
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

# Big export for plotting elsewhere...
saveRDS(list(OG_pleiotropy_Z = OG_pleiotropy_Z,
             tau_demo = tau_demo,
             centrality_combined = centrality_combined,
             pleiotropy_results_fig = pleiotropy_results_fig,
             climate_specific_Z = climate_specific_Z),
        paste0('outputs/',output_name,'_all_pleiotropy_figs_and_OG_Zscores.rds'))