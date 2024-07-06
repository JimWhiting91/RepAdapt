# Assess potential for sampling bias among datasets
# Script brings together LOO CV results with features of 
# datasets that might affect power
lib = c('data.table',
        'ggplot2',
        'dplyr',
        'tidyr',
        'vegan')
sapply(lib,library,character.only = T)

# Set some vars
run_name = "230321"
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"

# What species map to what genome?
species_genome_map = data.table(species =  c("Amaranthus tuberculatus", 
                                             "Arabidopsis halleri",     
                                             "Arabidopsis lyrata",      
                                             "Arabidopsis thaliana",    
                                             "Arabis alpina",           
                                             "Boechera stricta",        
                                             "Capsella rubella",        
                                             "Cardamine resedifolia", 
                                             "Eucalyptus albens",       
                                             "Eucalyptus magnificata",  
                                             "Eucalyptus sideroxylon",
                                             "Helianthus annuus",       
                                             "Helianthus argophyllus",  
                                             "Helianthus petiolaris",   
                                             "Medicago truncatula",     
                                             "Panicum hallii",          
                                             "Picea abies",             
                                             "Picea glaucaxengelmannii",
                                             "Picea obovata",           
                                             "Pinus contorta",          
                                             "Pinus sylvestris",        
                                             "Populus deltoides",       
                                             "Populus tremula",       
                                             "Populus trichocarpa",
                                             "Quercus petraea"),
                                genome = c('Atubercatus',
                                           'Ahalleri',
                                           'Alyrata',
                                           'Athaliana',
                                           'Aalpina',
                                           'Bstricta',
                                           'Crubella',
                                           'Ahalleri',
                                           'Egrandis',
                                           'Egrandis',
                                           'Egrandis',
                                           'Hannuus',
                                           'Hannuus',
                                           'Hannuus',
                                           'Mtruncatula',
                                           'Phallii',
                                           'Pabies',
                                           'Pabies',
                                           'Pabies',
                                           'Ptaeda',
                                           'Ptaeda',
                                           'Pdeltoides',
                                           'Ptremula',
                                           'Ptrichocarpa',
                                           'Qpetraea'))

# Fetch the cross-validation results
cv_res = readRDS(paste0("outputs/",output_name,"_picmin_results_SpeciesCrossValidation.rds"))
# Fetch the niche-breadth results
nb_res = fread(paste0('outputs/',output_name,'_gbif_niche_overlap.csv'))
# Fetch the original picmin results
true_picmin_outputs = readRDS(paste0("outputs/",output_name,"_picmin_results_doubleDS_pvals.rds"))
picmin_fdr = rbindlist(lapply(true_picmin_outputs,'[[',1))
picmin_fdr_signif = picmin_fdr[picmin_fdr < 0.5,]
picmin_tested_pvals = rbindlist(lapply(true_picmin_outputs,'[[',2))
# Fetch the Athal map
OG_map_Athal = data.table(readRDS(paste0("data/OG_map_Athal_",output_name,".rds")))

# Analyse the cross-validation results relative to the originals ----------

# Print the different number of FDR<0.5 results
for(i in 1:length(cv_res)){
  print(paste0(cv_res[[i]]$focal_species[1],': ',nrow(cv_res[[i]][picmin_fdr < 0.5,])))
}

# # What's the overlap of the retained OGs?
# species_retained_signif = lapply(cv_res,function(x) {
#   
#   # What OG were retained?
#   retained_OG = apply(x[picmin_fdr < 0.5,.(Orthogroup,climate_var)],1,paste,collapse = '-') |>
#     unique()
#   retained_OG
# })
# names(species_retained_signif) <- sapply(cv_res,function(x) x$focal_species[1])
# 
# sum(species_retained_signif$`Eucalyptus albens` %in% species_retained_signif$`Eucalyptus sideroxylon`)
# sum(species_retained_signif$`Eucalyptus sideroxylon` %in% species_retained_signif$`Eucalyptus albens`)

# What's the overlap of the removed OG from each set?
species_signif = lapply(cv_res,function(x) {
  
  # What OG were retained?
  retained_OG = apply(x[picmin_fdr < 0.5,.(Orthogroup,climate_var)],1,paste,collapse = '-') |>
    unique()
  # What OG have been removed from full dataset?
  original_OG = paste0(picmin_fdr_signif$Orthogroup,'-',picmin_fdr_signif$climate_var)
  # combined_OG = c(original_OG,retained_OG)
  changed_OG = original_OG[!original_OG %in% retained_OG]
  changed_OG
})
names(species_signif) <- sapply(cv_res,function(x) x$focal_species[1])

# make a matrix to show
overlap_mat = matrix(nrow = length(cv_res),ncol = length(cv_res))
for(i in 1:nrow(overlap_mat)){
  for(j in 1:ncol(overlap_mat)){
    if(i != j){
      overlap_mat[i,j] <- length(Reduce(intersect,species_signif[c(i,j)])) / length(species_signif[[i]])
      overlap_mat[j,i] <- length(Reduce(intersect,species_signif[c(i,j)])) / length(species_signif[[j]])
    }
  }
}
colnames(overlap_mat) <- rownames(overlap_mat) <- sapply(cv_res,function(x) x$focal_species[1])

cv_pairwise_heatmap = pheatmap::pheatmap(overlap_mat,
                                         main = 'Pairwise comparison of change in FDR<0.5 Orthogroups')
cv_pairwise_heatmap
# Plot out stacked bars of fdr signif hits
# Fetch a list of the number of FDR hits at different significance thresholds...
cross_val_signif_hits = lapply(cv_res,function(x){
  
  hits = x[picmin_fdr < 0.5,]
  hits$level = '<0.5'
  hits$level[hits$picmin_fdr < 0.4] = '<0.4'
  hits$level[hits$picmin_fdr < 0.3] = '<0.3'
  hits$level[hits$picmin_fdr < 0.2] = '<0.2'
  hits$level[hits$picmin_fdr < 0.1] = '<0.1'
  
  data.table(table(hits$level)) |>
    dplyr::rename('fdr_level' = 'V1') |>
    mutate(focal_species = x$focal_species[1])
  
}) |>
  rbindlist()

# Add to this the original hit numbers
true_hits = picmin_fdr_signif
true_hits$level = '<0.5'
true_hits$level[true_hits$picmin_fdr < 0.4] = '<0.4'
true_hits$level[true_hits$picmin_fdr < 0.3] = '<0.3'
true_hits$level[true_hits$picmin_fdr < 0.2] = '<0.2'
true_hits$level[true_hits$picmin_fdr < 0.1] = '<0.1'

true_signif_hits = data.table(table(true_hits$level)) |>
  dplyr::rename('fdr_level' = 'V1') |>
  mutate(focal_species = 'All Species')

# Combine and plot
fdr_05_cv = rbind(cross_val_signif_hits,
                  true_signif_hits) |>
  ggplot(aes(y = focal_species,x = N,fill = fdr_level)) +
  geom_bar(stat = 'identity',colour = 'black') +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = 'bottom') +
  scale_fill_manual(values = RColorBrewer::brewer.pal(5,'YlGnBu')) +
  geom_vline(xintercept = sum(true_signif_hits$N),linetype = 'dashed') +
  ggtitle('Cross-Validation (LOO) of\nSpecies on FDR <0.5 PicMin Results') +
  labs(x = 'N Signif PicMin Orthogroups (FDR <0.5)',
       y = 'Species Removed (LOO)',
       fill = 'FDR Level')

# Combine and plot at FDR <0.3
fdr_03_cv = rbind(cross_val_signif_hits,
                  true_signif_hits) |>
  subset(fdr_level %in% c('<0.1','<0.2','<0.3')) |>
  ggplot(aes(y = focal_species,x = N,fill = fdr_level)) +
  geom_bar(stat = 'identity',colour = 'black') +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = 'bottom') +
  scale_fill_manual(values = RColorBrewer::brewer.pal(5,'YlGnBu')) +
  geom_vline(xintercept = sum(true_signif_hits[fdr_level %in% c('<0.1','<0.2','<0.3'),N]),linetype = 'dashed') +
  ggtitle('Cross-Validation (LOO) of\nSpecies on FDR <0.3 PicMin Results') +
  labs(x = 'N Signif PicMin Orthogroups (FDR <0.3)',
       y = 'Species Removed (LOO)',
       fill = 'FDR Level')

# Combine these
fdr_bars = cowplot::plot_grid(fdr_05_cv,
                   fdr_03_cv,
                   ncol = 2,axis = 'tblr',align = 'h',
                   labels = c('A','B'),label_size = 32)

# Calculate the gbif distance ratio ---------------------------------------
# These are in km
gbif_dists = fread(paste0('outputs/',output_name,'_gbif_avg_distance.csv'))

# Calculate the distance ratio relative to the gbif data...
gbif_dists$distance_ratio = gbif_dists$avg_sampled_dist / gbif_dists$avg_gbif_distance

# Get some stats from original pvals --------------------------------------
# Pval files
pvals_file = paste0("outputs/GEA_res/run",run_name,"_",output_name,"_WZA_OG_PerGene_pvals.rds")
OG_pergene_pvals = readRDS(pvals_file)

# Orthogroup files
OG_map = readRDS(paste0("data/OF_OG_gea_gene_map_",output_name,".rds"))

# What were the tested orthogroup
tested_OG = unique(rbindlist(lapply(true_picmin_outputs,'[[',1))$Orthogroup)

# Loop through each species to find % of genome that was tested
species_genome_rep = lapply(1:nrow(species_genome_map),function(i){
  
  # Get codes of tested genes
  tested_genes = OG_pergene_pvals$mean_temp[species == species_genome_map$species[i] & Orthogroup %in% tested_OG,full_OF]
  # Get codes of all genes
  all_genes = OG_map[genome == species_genome_map$genome[i],full_OF]
  
  # Return the %
  data.table(species = species_genome_map$species[i],
             tested_percent = sum(all_genes %in% tested_genes)/length(all_genes))
}) |>
  rbindlist()


# Fetch the number of SNPs tested -----------------------------------------
# All this information is kept in Table S1
dataset_info = readxl::read_xlsx('tables/TableS1_dataset_info.xlsx',skip = 1) |>
  as.data.table()
dataset_info$Species = gsub('Picea glauca x engelmannii',
                            'Picea glaucaxengelmannii',
                            dataset_info$Species)

# Average the values for when multiple datasets were combined to one
dataset_info = dataset_info[,.(`SNP N` = mean(.SD$`SNP N`),
                               `Ind N` = mean(.SD$`Ind N`),
                               `Sampling Site N` = mean(.SD$`Sampling Site N`)),by = Species] |>
  mutate(`Ind:Site Ratio` = `Ind N`/`Sampling Site N`)

# Build a combined data table for a heatmap -------------------------------
# Results from CV LOO
all_res = lapply(cv_res,function(x){
  out = data.table(cv_signif = nrow(x[picmin_fdr < 0.5,]),
                   species = x$focal_species[1],
                   cv_change = nrow(x[picmin_fdr < 0.5,]) - nrow(picmin_fdr_signif))
}) |>
  rbindlist() |>
  dplyr::rename(Species = species)

# Add the niche breadth
all_res = merge(all_res,nb_res,
                by.x = 'Species',by.y = 'species')
# Add the geographic distance
all_res = merge(all_res,gbif_dists[,.(species,distance_ratio)],
                by.x = 'Species',by.y = 'species')
# Add the % of genes tested
all_res = merge(all_res,species_genome_rep,
                by.x = 'Species',by.y = 'species')
# Add the number of SNPs
# Add the number of individuals
# Add the number of locations
all_res = merge(all_res,dataset_info) |>
  dplyr::rename(`Global Climate Niche Overlap` = wmean_overlap,
                `Global Range Distance Ratio` = distance_ratio,
                `% of Genes Tested` = tested_percent)

# Calculate the nonparametric correlation between each feature and the LOO res
loo_corrs = cor(all_res[,3:ncol(all_res)],method = 'spearman')
loo_corrs = reshape2::melt(loo_corrs['cv_change',colnames(loo_corrs) != 'cv_change'])
loo_corrs$Feature = rownames(loo_corrs)

# What are the signifs of the corrs?
feature_cors = lapply(4:ncol(all_res),function(i){
  cor.test(all_res$cv_change,
           all_res[[i]],
           method = 'spearman')
})

# Plot them
loo_corr_plot = ggplot(loo_corrs,aes(y = Feature,x = value,fill = value)) +
  geom_bar(stat = 'identity',colour = 'black') +
  scale_fill_gradient2() +
  theme_minimal() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        legend.position = 'none') +
  labs(x = expression(Spearman~rho)) +
  ggtitle('Associations between dataset features\nand LOO CV results')
loo_corr_plot

# Make a heatmap just showing the CV change alongside the features
# Scale these to some common scale
scale_all_res = all_res
for(i in 3:ncol(scale_all_res)){
  scale_all_res[[i]] <- scale(scale_all_res[[i]])
}
species_feature_heatmap = scale_all_res |>
  dplyr::select(-cv_signif) |>
  dplyr::rename(`LOO CV Change` = cv_change) |>
  mutate(species_F = factor(Species,levels = Species[order(`LOO CV Change`)])) |>
  reshape2::melt() |>
  mutate(facet_var = ifelse(variable == 'LOO CV Change',
                            'Result','Dataset Feature')) |>
  ggplot(aes(y = species_F,x = variable,fill = value)) +
  geom_tile() +
  scale_fill_gradient2() +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12,angle = 45,hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_blank(),
        legend.position = 'top') +
  labs(y = 'Species',
       fill = 'Scaled Value')

species_feature_fig = cowplot::plot_grid(species_feature_heatmap,
                   loo_corr_plot,
                   ncol = 2,axis = 'tblr',align = 'h',
                   labels = c('C','D'),
                   label_size = 32)


# Extended Data plot ------------------------------------------------------
pdf('figs/FigureSX_LOO_res.pdf',width = 13,height = 14)
cowplot::plot_grid(fdr_bars,
                   species_feature_fig,
                   nrow = 2)
dev.off()


# Redundancy Analysis of results ------------------------------------------
# Do an RDA to model collectively the seven features vs the cv change 
y_var = all_res$cv_change
x_var = data.frame(scale(all_res[,4:ncol(all_res)]))
rownames(x_var) = all_res$Species
# Build the model
loo_rda = rda(y_var ~ .,data = x_var)

summary(loo_rda)

