####  Analysis of associations between duplication in gene trees and picmin results
lib = c('ghibli','ggtree','adephylo','ape','dplyr',"cowplot","ggridges","viridis","tidyr","ape","DHARMa","sjPlot","qvalue","ggplot2","pbmcapply","biomaRt","org.At.tair.db","ExpressionAtlas","S4Vectors", "IRanges", "GenomicRanges", "SummarizedExperiment","readr","data.table")
sapply(lib,library,character.only=T)
n_cores = 6
source("R/repadapt_functions.R")

# Functions ---------------------------------------------------------------
prune_duplicates = function(dups_file = OF_duplicates_HQ,
                            OG,
                            OG_pvals = OG_pvals,
                            species_nodes = species_nodes){
  
  # Which genomes did we test...
  species_tested = unique(OG_pvals[Orthogroup == OG,species])
  genomes_tested = species_genome_map[species %in% species_tested,genome]
  
  # Get rid of cases where the genome or node isn't tested
  nodes_to_keep = unique(species_nodes[species %in% genomes_tested,node])
  
  # Return
  OF_duplicates_HQ[Orthogroup == OG & `Species Tree Node` %in% c(genomes_tested,nodes_to_keep),]
}


# Prepare metadata --------------------------------------------------------

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

# Duplication branches... -------------------------------------------------
OF_dir = "outputs/orthology/Results_221213_18_genomes_Ptaeda_isoforms_removed/"
OF_duplicates = fread(paste0(OF_dir,"/Gene_Duplication_Events/Duplications.tsv"))
OF_duplicates_HQ = OF_duplicates[Support > 0.7,]
OF_tree_dir = paste0(OF_dir,"/Resolved_Gene_Trees/")


# We want to filter the duplications for cases where we removed there being lots and lots of paralogs and only retain duplications with downstream genomes we looked at...
# Use the species tree to filter
OF_species_tree = read.tree(paste0(OF_dir,"/Species_Tree/SpeciesTree_rooted_node_labels.txt"))
all_nodes = listTips(OF_species_tree)
species_nodes = rbindlist(lapply(1:length(all_nodes),function(i){
  out = data.table(node = names(all_nodes)[i],
                   species = names(all_nodes[[i]]))
}))

OF_duplicates_HQ_pruned = rbindlist(pbmclapply(unique(OF_duplicates_HQ$Orthogroup),prune_duplicates,
                                               dups_file = OF_duplicates_HQ,
                                               OG_pvals = OG_pvals,
                                               species_nodes = species_nodes,
                                               mc.cores = n_cores))

# Plot an example duplication tree ----------------------------------------
# Just make up some tree with dups in it...
myTree <- ape::read.tree(text='(S1,(((S4,(S3,S3)D2),S2),((S3,S4),(S2,S2)D3))D1);')
myTree$node.label = c("","italic('D1')",'','',"italic('D2')",'','',"italic('D3')")
branches <- list(S1=1,S2=c(5,8,9),S3=c(3,4,6),S4=c(2,7))
myTree <- groupOTU(myTree,branches)
dup_tree_example = ggtree(myTree) +
  # geom_text(aes(label=node)) +
  geom_nodelab(hjust = -0.1,size= 6,parse = T) +
  geom_tiplab(aes(color=group),size = 10,)+
  scale_color_manual(values=c(S1 = "#26432FFF",S2= "#4D6D93FF", S3 = "#6FB382FF", S4 = "#E48C2AFF")) +
  theme(legend.position = 'none') + 
  ggplot2::xlim(0,6)
dup_tree_example

# The number of duplications per OG ---------------------------------------
# All duplications
OG_dupN = data.table(table(OF_duplicates_HQ_pruned$Orthogroup))
colnames(OG_dupN)[1] = "Orthogroup"
OG_dupN = OG_dupN[Orthogroup %in% picmin_fdr$Orthogroup,]

# Add on zeroes...
OG_dupN = rbind(OG_dupN,
                data.table(Orthogroup = unique(picmin_fdr[!Orthogroup %in% OG_dupN$Orthogroup,Orthogroup]),
                           N = 0))

# All species-specific dups...
OG_dupN_species = data.table(table(OF_duplicates_HQ_pruned[`Species Tree Node` %in% unique(species_nodes$species),Orthogroup]))
colnames(OG_dupN_species)[1] = "Orthogroup"
OG_dupN_species = OG_dupN_species[Orthogroup %in% picmin_fdr$Orthogroup,]

OG_dupN_species = rbind(OG_dupN_species,
                        data.table(Orthogroup = unique(picmin_fdr[!Orthogroup %in% OG_dupN_species$Orthogroup,Orthogroup]),
                                   N = 0))

# Also count up the number of single copy species...
singlecopy_species = OG_pvals[Orthogroup %in% names(OG_counts)[OG_counts == max(OG_counts)] & climate_var == 'mean_temp',
                              .(single_copy_N = sum(.SD$Ngenes_per_species == 1)),by = Orthogroup]

# Merge with picmin minp
dup_merge = merge(OG_dupN,OG_dupN_species,by = 'Orthogroup')
colnames(dup_merge) = c("Orthogroup","all_dups_N","species_dups_N")
dup_merge = merge(dup_merge,singlecopy_species)
dup_merge = merge(dup_merge,picmin_deciles)

# Plot these...
plot_dups = melt(dup_merge)[variable != 'PicMin_minP',.(decile_mean = mean(value)),by = .(picmin_decile,variable)][order(picmin_decile)]
plot_dups$variable_labs = stringr::str_to_title(gsub("_"," ",plot_dups$variable))
# Also add means per variable...
means_to_plot = data.table(variable_labs = c("All Dups N","Species Dups N","Single Copy N"),
                           mean_vals = c(mean(dup_merge$all_dups_N),
                                         mean(dup_merge$species_dups_N),
                                         mean(dup_merge$single_copy_N)))

# Repeat the analysis above but using picmin deciles from random results --------
#####
# NOTE: This analysis takes LOTS of memory, only run if necessary
run_randomised_dups_analysis = FALSE
#####
if(run_randomised_dups_analysis){
  random_picmin_res = readRDS(paste0("outputs/",output_name,"_randomised_pvals_picmin_res.rds"))
  random_picmin_res = data.table(separate(random_picmin_res,col = 'test_group', into = c("Orthogroup","climate_var"),sep = ':'))
  
  # Repeat the decile approach
  random_picmin_deciles = random_picmin_res[,.(random_minP = min(picmin_p_adj)),by = .(Orthogroup,iteration)]
  random_picmin_deciles = rbindlist(lapply(unique(random_picmin_deciles$iteration),function(iter){
    tmp = random_picmin_deciles[iteration == iter,]
    tmp$random_picmin_decile = cut_number(tmp$random_minP,n = 10)
    tmp$random_picmin_decile_order = 0
    for(f in 1:length(levels(tmp$random_picmin_decile))){
      tmp[random_picmin_decile == levels(tmp$random_picmin_decile)[f],random_picmin_decile_order := f]
    }
    
    # Merge and calc stats here...
    tmp_merge = merge(dup_merge,tmp)
    # Within each decile set, summarise
    out = tmp_merge[,.(all_dups_N = mean(all_dups_N),
                       species_dups_N = mean(species_dups_N),
                       single_copy_N = mean(single_copy_N)),by = random_picmin_decile_order]
    out$iteration = iter
    out
  }))
  
  # For each decile group we can now take means and quantiles
  random_picmin_deciles$random_picmin_decile_order = factor(random_picmin_deciles$random_picmin_decile_order)
  random_to_plot = melt(random_picmin_deciles)[variable != "iteration",][,.(random_mean = mean(value),
                                                                            random_q5 = quantile(value,probs = 0.05),
                                                                            random_q95 = quantile(value,probs = 0.95),
                                                                            min_val = min(value),
                                                                            max_val = max(value)),by = .(variable,random_picmin_decile_order)]
  
  # Add on to these the original picmin deciles...
  picmin_decile_groups = sort(unique(plot_dups$picmin_decile))
  random_to_plot = random_to_plot[order(random_picmin_decile_order),] |>
    mutate(picmin_decile = rep(picmin_decile_groups,each = 3))
  random_to_plot$variable_labs = stringr::str_to_title(gsub("_"," ",random_to_plot$variable))
  
  # Plot all of the duplication results -------------------------------------
  duplications_results_fig = ggplot(plot_dups,aes(y = picmin_decile,x = decile_mean)) +
    geom_point(data = random_to_plot,aes(y = picmin_decile,x = random_mean),colour = 'black',size = 4) +
    geom_segment(data = random_to_plot,aes(y = picmin_decile,yend = picmin_decile,
                                           x = min_val,xend = max_val)) +
    geom_point(size = 4,colour = 'red3',shape = 17) +
    facet_wrap(~variable_labs,nrow = 1,scale = 'free_x') +
    scale_y_discrete(limits = rev) +
    labs(y = "Evidence for repeatability\n(min PicMin p-value)",x = "Mean duplication index per-decile") +
    theme_bw() +
    theme(legend.position = 'right',
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12),
          strip.text = element_text(size = 12),
          strip.background = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 12))
  duplications_results_fig
}

# Duplications per species ------------------------------------------------
# Here we want to separate repeatable orthogroups into 'contributing' and not, and compare the distribution of duplications?
convergent_picmin_res = picmin_fdr[picmin_fdr < 0.5,]
contributing_convergent_OG_pvals = rbindlist(lapply(1:nrow(convergent_picmin_res),function(x){
  
  # What's the cutoff
  cutoff_tmp = convergent_picmin_res$config_est[x]
  tmp = OG_pvals[climate_var == convergent_picmin_res$climate_var[x] &
                   Orthogroup == convergent_picmin_res$Orthogroup[x],][order(epval_final),]
  p_cutoff = tmp$epval_final[cutoff_tmp]
  if(p_cutoff > 0.1){
    p_cutoff = 0.1
  }
  tmp = tmp[epval_final < p_cutoff,]
  tmp$convergent_index = x
  tmp
}))

# What's the average paralog N in these
contributing_mean_dups = contributing_convergent_OG_pvals[,.(mean_dups = mean(Ngenes_per_species),
                                                             prop_dups = 1 - (sum(.SD$Ngenes_per_species == 1)/nrow(.SD))),by = convergent_index]

# Exclude those cases and repeat 10,000 times...
nonconvergent_picmin_res = picmin_fdr[picmin_fdr > 0.5,]
OG_pvals$clim_OG = paste0(OG_pvals$Orthogroup,":",OG_pvals$climate_var)
random_contributing_perms = rbindlist(pbmclapply(1:10000,function(iter){
  set.seed(iter)
  random_res = nonconvergent_picmin_res[sample(1:nrow(nonconvergent_picmin_res),nrow(convergent_picmin_res)),]
  climOG_tmp = paste0(random_res$Orthogroup,":",random_res$climate_var)
  tmp_pvals = OG_pvals[clim_OG %in% climOG_tmp & epval_final < 0.1,]
  
  # Calculate per OG stats...
  tmp_pvals_stats = tmp_pvals[,.(mean_dups = mean(Ngenes_per_species),
                                 prop_dups = 1 - (sum(.SD$Ngenes_per_species == 1)/nrow(.SD))),by = clim_OG]
  out = data.table(iter = iter,
                   mean_dups_per_contributing = mean(tmp_pvals_stats$mean_dups),
                   mean_props_duplicated = mean(tmp_pvals_stats$prop_dups))
},mc.cores = n_cores))

# Calculate empirical p-val
# Total number of dups
contributing_dups_p = qvalue::empPvals(mean(contributing_mean_dups$mean_dups),random_contributing_perms$mean_dups_per_contributing)
# Proportion duplicated
contributing_prop_dupped_p = qvalue::empPvals(mean(contributing_mean_dups$prop_dups),random_contributing_perms$mean_props_duplicated)

# Make figs of these
plot_contributing = data.table(random_perms = c(random_contributing_perms$mean_dups_per_contributing,
                                                random_contributing_perms$mean_props_duplicated),
                               variable = rep(c("Mean Duplication N","Mean Proportion of\nDuplicated Genes"),each = nrow(random_contributing_perms)))
obs_to_plot = data.table(obs = c(mean(contributing_mean_dups$mean_dups),mean(contributing_mean_dups$prop_dups)),
                         variable = c("Mean Duplication N","Mean Proportion of\nDuplicated Genes"),
                         signif = ifelse(c(contributing_dups_p,contributing_prop_dupped_p) < 0.05,"Signif","No"))

contributing_species_fig = ggplot(plot_contributing,aes(x = random_perms)) +
  geom_histogram(colour = 'gray50') +
  facet_wrap(~variable,nrow = 1,scales = 'free') +
  theme_bw() +
  labs(y = 'Count',x = 'Random Permuted Mean') +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 13),
        strip.text = element_text(size = 13),
        legend.position = 'none',
        strip.background = element_blank(),
        title = element_text(size = 14)) +
  geom_vline(data = obs_to_plot,aes(xintercept = obs,colour = signif),size = 2) +
  scale_colour_manual(values = c("Signif" = "green4","No" = "black")) +
  ggtitle("Contributing species within RAOs")
contributing_species_fig


# Is duplication linked to pleiotropy? ------------------------------------
# Fetch all of the per Orthogroup pleiotropy metrics...
pleiotropy_outputs = readRDS(paste0('outputs/',output_name,'_all_pleiotropy_figs_and_OG_Zscores.rds'))

# Merge the pleiotropy Z scores with OG Z scores...
all_OG_Z = merge(OG_dupN,pleiotropy_outputs$OG_pleiotropy_Z,by = 'Orthogroup')

# We also want to make some random nonsense variable that shares the orthogroup structure to confirm that this isn't related...
# Use the tau data structure...
Athal_tau_res = readRDS("outputs/Athal_tau_data.rds")
OG_tau = merge(Athal_tau_res,OG_map_Athal[,.(TAIR_gene,Orthogroup)])
# Loop over 100 runs
random_corrs = rbindlist(lapply(1:100,function(x){
  print(x)
  set.seed(x)
  
  # This is just an empirical nonsense set of random_pvals
  random_pvals = runif(nrow(OG_tau))
  OG_tau$tau_p = empPvals(-random_pvals,-random_pvals)
  random_OG_Z = condenseToOGPvals(na.omit(OG_tau),"Orthogroup","tau_p")
  random_OG_Z$rand_Z = qnorm(1 - random_OG_Z$minP_DS_epvalue)
  
  # Merge this and get the corrs...
  random_OG_Z_merge = merge(all_OG_Z,random_OG_Z[,.(Orthogroup,rand_Z)])
  random_corrs = cor(random_OG_Z_merge[,2:ncol(random_OG_Z_merge)],method = 'spearman')[,"rand_Z"]
  out = data.table(variable = names(random_corrs),
                   random_corrs)
}))
random_corrs_means = random_corrs[,.(mean_corr = mean(random_corrs)),by = variable]

# Add these to the corr mat
tmp_mat = cor(all_OG_Z[,2:ncol(all_OG_Z)])
colnames(tmp_mat)[1] <- rownames(tmp_mat)[1] <- "Duplication N"
functional_corr_mat = matrix(ncol = ncol(tmp_mat) + 1,
                             nrow = nrow(tmp_mat) + 1)
functional_corr_mat[1:ncol(tmp_mat),1:ncol(tmp_mat)] = tmp_mat
functional_corr_mat[,ncol(functional_corr_mat)] = random_corrs_means$mean_corr
functional_corr_mat[nrow(functional_corr_mat),] = random_corrs_means$mean_corr
colnames(functional_corr_mat) <- rownames(functional_corr_mat) <- stringr::str_to_title(gsub("_"," ",c(colnames(tmp_mat),"Random")))
colnames(functional_corr_mat) <- rownames(functional_corr_mat) <- gsub("Node ","",rownames(functional_corr_mat))
colnames(functional_corr_mat) <- rownames(functional_corr_mat) <- gsub("Tau Z","Exp Breadth Z",rownames(functional_corr_mat))

# Melt this and plot...
functional_corr_dd = data.table(melt(functional_corr_mat))
functional_corr_dd$Var1_F = factor(functional_corr_dd$Var1,levels=colnames(functional_corr_mat))
functional_corr_dd$Var2_F = factor(functional_corr_dd$Var2,levels=colnames(functional_corr_mat))
functional_corrs = ggplot(functional_corr_dd,aes(Var1_F,Var2_F,fill = value)) +
  geom_tile(colour = 'black') +
  scale_fill_gradient2(low = 'blue4',mid = 'white',high = 'red4',midpoint = 0,
                       limits = c(-1,1),breaks = c(-1,0,1)) +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1),
        axis.title = element_blank(),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = 'top') +
  labs(fill = expression(Spearman~rho))
functional_corrs

# Fetch some specific correlations
cor.test(all_OG_Z$At_node_degree_Z,all_OG_Z$N,method = 'spearman')
cor.test(all_OG_Z$tau_Z,all_OG_Z$N,method = 'spearman')

# Plot the pleiotropy results and duplication results separately...
# Dup results
dups_final_fig = cowplot::plot_grid(dup_tree_example,
                                    duplications_results_fig,
                                    contributing_species_fig,
                                    ncol = 3,rel_widths = c(1,2,1.5),
                                    labels = "AUTO",label_size = 32)
pdf('figs/FigureSX_all_duplication_results.pdf',width = 15,height = 4)
dups_final_fig
dev.off()

# Pleiotropy...
functional_patch1 = cowplot::plot_grid(pleiotropy_outputs$tau_demo,
                                       pleiotropy_outputs$centrality_combined +
                                         theme_minimal(),
                                       nrow = 1,labels = c("A","C"),label_size = 32,hjust = c(0,1),rel_widths = c(0.6,1))
functional_patch1

functional_patch2 = cowplot::plot_grid(pleiotropy_outputs$pleiotropy_results_fig,
                                       cowplot::plot_grid(pleiotropy_outputs$climate_specific_Z,
                                                          functional_corrs,
                                                          ncol = 2,labels = c("D","E"),label_size = 32,rel_widths = c(1.5,1)),
                                       ncol = 1,labels = c("B",""),label_size = 32,rel_heights = c(1,1),vjust = c(1.4,0))
functional_patch2

# Merge
pdf('figs/Figure5_pleiotropy_and_convergence_results.pdf',width = 14,height = 12)
cowplot::plot_grid(functional_patch1,
                   functional_patch2,
                   ncol = 1,
                   rel_heights = c(1,2.8))
dev.off()
