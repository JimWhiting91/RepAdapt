# Test script for exploring STRINGdb analyses...
lib <- c("ggdendro","GOSemSim","dplyr","simplifyEnrichment","igraph","ggridges","ggrepel","cowplot","ggplot2","data.table","tidyr","pbmcapply","ggtree","ape","viridis","STRINGdb","qvalue","dplyr","biomaRt")
sapply(lib,library,character.only=T)
source("R/repadapt_functions.R")

# Which run name are we looking at?
run_name = "230321"
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"
OG_map_Athal = data.table(readRDS(paste0("data/OG_map_Athal_",output_name,".rds")))
OG_map = data.table(readRDS(paste0("data/OF_OG_gea_gene_map_",output_name,".rds")))
n_cores = 6

# Functions ---------------------------------------------------------------
boot_OG_interactions = function(string_vector,quality){
  # Set cutoff
  if(quality == 'low'){
    cutoff = 150
  } else if(quality == 'med'){
    cutoff = 400
  } else if(quality == 'high'){
    cutoff = 700
  } else if(quality == 'highest'){
    cutoff = 900
  } else {
    cutoff = 0
  }
  # Fetch the interactions of all genes
  tmp_interactions = data.table(Athal_string_db$get_interactions(string_vector))
  nrow(unique(tmp_interactions[combined_score > cutoff,]))
}

climate_var_title = function(x){
  stringr:::str_to_title(gsub("_"," ",x))
}

# Fetch inputs ----------------------------------
picmin_outputs = readRDS(paste0("outputs/",output_name,"_picmin_results_doubleDS_pvals.rds"))
picmin_fdr = rbindlist(lapply(picmin_outputs,'[[',1))
OG_pvals = rbindlist(lapply(picmin_outputs,'[[',2))

# Prep STRING -------------------------------------------------------------
ncbi_species <- read.csv("metadata/filtered_species_list_ncbi.txt",header=F)
colnames(ncbi_species) <- c("species","ncbi")

# For now, just look at Athal...
Athal_string_db <- STRINGdb$new(version="11.5", species=ncbi_species[grep("thaliana",ncbi_species$species),"ncbi"],score_threshold=200, input_directory="")

# First fetch the all of the string aliases...
string_aliases = data.table(Athal_string_db$get_aliases(takeFirst = F))

# For the orthogroups we have studied with picmin, fetch their aliases...
to_fetch = unique(picmin_fdr$Orthogroup)

# For the first of these, we can simply assume their string codes...
OG_string_codes = OG_map_Athal[Orthogroup %in% to_fetch,.(Orthogroup,TAIR_gene),]
OG_string_codes = OG_string_codes[grep("AT",OG_string_codes$TAIR_gene),]
OG_string_codes$STRING_id = paste0(ncbi_species[grep("thaliana",ncbi_species$species),"ncbi"],".",OG_string_codes$TAIR_gene,".1")

# Resolve the ones that don't match...
to_resolve = OG_string_codes$STRING_id[!OG_string_codes$STRING_id %in% string_aliases$STRING_id]
true_strings = unique(string_aliases$STRING_id)
OG_string_codes_fixed = pbmclapply(to_resolve,function(tmp_string){
  
  new_string = grep(OG_string_codes[STRING_id == tmp_string,TAIR_gene],true_strings,value = T)
  if(length(new_string) != 1){
    return(NULL)
  } else {
    out = OG_string_codes[STRING_id == tmp_string,]
    out$STRING_id = new_string
    return(out)
  }
},mc.cores = n_cores)

# FINALLY merge the good strings with the resolved strings...
OG_string_codes = rbind(OG_string_codes[OG_string_codes$STRING_id %in% string_aliases$STRING_id,],
                        rbindlist(OG_string_codes_fixed))

# Are climate convergent OG more connected than neutral?-------------------
# First we want to do 1000 bootstraps where we retain a single gene per orthogroup and assess connectivity of the network...
bootN = 1000
fdr_thresh = 0.5
string_quality = 'med'

# We'll test all OG together, then separate into climate and precipitation
climate_var_list = list(all_vars = unique(picmin_fdr$climate_var),
                        temp_vars = grep("temp|tmax|iso|diurnal",unique(picmin_fdr$climate_var),value = T),
                        precip_vars = grep("temp|tmax|iso|diurnal",unique(picmin_fdr$climate_var),value = T,invert = T))

climate_OG_convergent = list(all_vars = unique(picmin_fdr[picmin_fdr < fdr_thresh,Orthogroup]),
                             temp_vars = unique(picmin_fdr[picmin_fdr < fdr_thresh & climate_var %in% climate_var_list$temp_vars,Orthogroup]),
                             precip_vars = unique(picmin_fdr[picmin_fdr < fdr_thresh & climate_var %in% climate_var_list$precip_vars,Orthogroup]))

climate_bootstrap_connect = lapply(names(climate_OG_convergent),function(var){
  message(paste0(">>> STARTING: ",var))
  
  boot_genes_to_test = rbindlist(pbmclapply(1:bootN,function(x){
    set.seed(x)
    # Draw our random genes from orthogroups...
    random_OG_genes = na.omit(OG_string_codes)[Orthogroup %in% climate_OG_convergent[[var]],][,.(STRING_id = sample(.SD$STRING_id,1)), by = Orthogroup]
    out = data.table(random_OG_genes,
                     iter = x)
  },mc.cores = n_cores))
  
  # Now identify the number of interactions between all string IDs in the boots...
  message(paste0(">>> Counting interactions..."))
  start = Sys.time()
  bootstrap_out = boot_genes_to_test[,.(interactionN = boot_OG_interactions(.SD$STRING_id,quality = string_quality)),by = iter]
  message(paste0(">>> Finished in ",round(Sys.time() - start)," seconds"))
  bootstrap_out$climate_var = var
  bootstrap_out
})

rbindlist(climate_bootstrap_connect)[,mean(interactionN),by = climate_var]


# Then repeat this but randomise it over randomN draws of Orthogroups...
randomN = 10000
if(!file.exists(paste0("outputs/",output_name,"_STRING_randomperms.rds"))){
  random_bootstrap_connect = lapply(names(climate_OG_convergent),function(var){
    message(paste0(">>> STARTING: ",var))
    
    random_genes_to_test = rbindlist(pbmclapply(1:randomN,function(x){
      set.seed(x)
      # Draw our random genes from orthogroups...
      true_OG = unique(na.omit(OG_string_codes)[Orthogroup %in% climate_OG_convergent[[var]],Orthogroup])
      random_OG = sample(unique(na.omit(OG_string_codes)[!Orthogroup %in% true_OG,Orthogroup]),length(true_OG))
      random_OG_genes = na.omit(OG_string_codes)[Orthogroup %in% random_OG,][,.(STRING_id = sample(.SD$STRING_id,1)), by = Orthogroup]
      out = data.table(random_OG_genes,
                       iter = x)
    },mc.cores = n_cores))
    
    # Now identify the number of interactions between all string IDs in the boots...
    message(paste0(">>> Counting interactions..."))
    start = Sys.time()
    bootstrap_out = random_genes_to_test[,.(interactionN = boot_OG_interactions(.SD$STRING_id,quality = string_quality)),by = iter]
    message(paste0(">>> Finished in ",round(Sys.time() - start)," seconds"))
    bootstrap_out$climate_var = var
    return(bootstrap_out)
  })
  saveRDS(random_bootstrap_connect,paste0("outputs/",output_name,"_STRING_randomperms.rds"))
} else {
  random_bootstrap_connect = readRDS(paste0("outputs/",output_name,"_STRING_randomperms.rds"))
}

# Get the means to plot and plot these relative to the whole distributions...
means_to_plot = rbindlist(climate_bootstrap_connect)[climate_var %in% names(climate_OG_convergent),
                                                     .(boot_mean = mean(interactionN)),by = climate_var]

# Get empPvals for these...
climate_interactions_epvals = rbindlist(lapply(names(climate_OG_convergent),function(var){
  epval =  empPvals(means_to_plot[climate_var == var,boot_mean],
                    rbindlist(random_bootstrap_connect)[climate_var == var,interactionN])
  out = data.table(climate_var = var,
                   interactionN_epval = epval)
}))
climate_interactions_epvals = climate_interactions_epvals[order(interactionN_epval),]
climate_interactions_epvals

random_to_plot = rbindlist(random_bootstrap_connect)
random_to_plot$climate_var = climate_var_title(random_to_plot$climate_var)
means_to_plot$climate_var = climate_var_title(means_to_plot$climate_var)

random_to_plot$climate_var_F = factor(random_to_plot$climate_var,levels =  climate_var_title(climate_interactions_epvals$climate_var))
means_to_plot$climate_var_F = factor(means_to_plot$climate_var,levels = climate_var_title(climate_interactions_epvals$climate_var))

# Mark signif
means_to_plot$boot_signif = "No"
means_to_plot[climate_var %in% climate_var_title(climate_interactions_epvals[interactionN_epval < 0.05,climate_var]),boot_signif := "Yes"]

# Also to this, add the number of orthogroups that we originally had to draw...
# convergent_OG_N = data.table(table(picmin_fdr[picmin_fdr < 0.5,climate_var]))
convergent_OG_N = data.table(climate_var = names(climate_OG_convergent),
                             N = sapply(climate_OG_convergent,length))
convergent_OG_N$climate_var = climate_var_title(convergent_OG_N$climate_var)
random_to_plot = merge(random_to_plot,convergent_OG_N[,.(N,climate_var)])

# Alternative figure for string perm res
random_to_plot2 = random_to_plot[,.(random_mean = mean(interactionN),
                                    q5 = quantile(.SD$interactionN,0.05),
                                    q95 = quantile(.SD$interactionN,0.95),
                                    N = unique(.SD$N)),by = climate_var_F]
signif_marks = means_to_plot[climate_var %in% climate_var_title(climate_interactions_epvals[interactionN_epval < 0.05,climate_var]),.(climate_var_F,boot_mean)]
signif_marks$x_val = max(signif_marks$boot_mean) + 0.5*max(signif_marks$boot_mean)
signif_marks$signif = "*"
for(i in 1:nrow(signif_marks)){
  pval_tmp = climate_interactions_epvals[climate_var_title(climate_var) == signif_marks$climate_var_F[i],interactionN_epval]
  if(pval_tmp < 0.001){
    signif_marks$signif[i] = "***"
  } else if(pval_tmp < 0.01){
    signif_marks$signif[i] = "**"
  }
}

# Make new titles
new_titles = c("All Variables","Temperature","Precipitation")
old_titles = c("All Vars","Temp Vars","Precip Vars")
for(i in 1:length(new_titles)){
  random_to_plot2$climate_var_F = gsub(old_titles[i],new_titles[i],random_to_plot2$climate_var_F)
  means_to_plot$climate_var_F = gsub(old_titles[i],new_titles[i],means_to_plot$climate_var_F)
  signif_marks$climate_var_F = gsub(old_titles[i],new_titles[i],signif_marks$climate_var_F)
}
random_to_plot2$climate_var_F = factor(random_to_plot2$climate_var_F,levels = new_titles)
means_to_plot$climate_var_F = factor(means_to_plot$climate_var_F,levels = new_titles)
signif_marks$climate_var_F = factor(signif_marks$climate_var_F,levels = new_titles)

string_perm_res_dots = ggplot(random_to_plot2,aes(x = random_mean,y = climate_var_F)) +
  geom_point(data = means_to_plot,aes(y = climate_var_F,x = boot_mean),shape = 17, colour = "forestgreen",size = 7) +
  # geom_point(data = means_to_plot[boot_signif == "No",],aes(y = climate_var_F,x = boot_mean),shape = 17, colour = "black",size = 4) +
  geom_text(data = signif_marks,aes(x = x_val,label = signif,y = climate_var_F),size = 10) +
  geom_segment(aes(x = q5,xend = q95,y = climate_var_F,yend = climate_var_F),size = 1.5) +
  geom_point(aes(colour=N), size=6) + 
  geom_point(shape = 1,size = 6,colour = "black")  +
  scale_y_discrete(limits=rev) +
  scale_colour_viridis(option = "C") +
  theme_minimal() +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        legend.position = "top",
        legend.title = element_text(size = 14),
        legend.text  = element_text(size = 12)) +
  xlim(0,1.25*(max(means_to_plot$boot_mean))) +
  labs(x = "Within-network interactions",colour = "Orthogroup N") +
  guides(colour = guide_colourbar(title.position="top", title.hjust = 0.5))
string_perm_res_dots


# Enrichment within networks ----------------------------------------------
ensembl <- useMart(biomart = "plants_mart",host="https://plants.ensembl.org")
athal_biomart <- useDataset("athaliana_eg_gene",mart=ensembl)
athal_universe <- getBM(attributes = c("tair_locus","ensembl_gene_id","chromosome_name","start_position","end_position",
                                       "go_id","name_1006"),
                        mart=athal_biomart)

# Merge with Orthogroups...
tmp = data.table(OG_map_Athal)
tmp$ensembl_gene_id = tmp$TAIR_gene
athal_universe_OG = data.table(merge(athal_universe,tmp[,.(ensembl_gene_id,Orthogroup)]))
colnames(athal_universe_OG)[which(colnames(athal_universe_OG) == "name_1006")] = "go_name"

# Set up unique OG-GO groups...
OG_GO = unique(athal_universe_OG[,.(go_id,go_name,Orthogroup)])
OG_GO = OG_GO[Orthogroup %in% unique(picmin_fdr$Orthogroup),]

# For each of our significant networks, do GO enrichment...
climate_enrichment = lapply(climate_OG_convergent,GO_enrichment_orthogroups,
                            GO_universe = OG_GO,
                            n_cores = 4)

# Plot these...
AtGO <- godata('org.At.tair.db', ont="BP")
var_titles = c("All","Temperature","Precipitation")
GO_species_heat_list = lapply(1:length(climate_enrichment),function(x){
  print(x)
  # Fetch the relevant OG
  go_tmp = climate_enrichment[[x]]$go_id
  clim_tmp = names(climate_enrichment)[x]
  OG_tmp = climate_OG_convergent[[clim_tmp]]
  
  # Match these with the go terms we are looking at
  OG_GO_tmp = OG_GO[Orthogroup %in% OG_tmp &
                      go_id %in% go_tmp,]
  
  # Fetch the relevant species now using the orthogroups
  OG_pvals_tmp = OG_pvals[climate_var %in% climate_var_list[[names(climate_OG_convergent)[x]]],]
  picmin_fdr_tmp = picmin_fdr[climate_var %in% climate_var_list[[names(climate_OG_convergent)[x]]] &
                                picmin_fdr < fdr_thresh &
                                Orthogroup %in% OG_GO_tmp$Orthogroup,]
  
  species_contributing_tmp = rbindlist(lapply(1:nrow(picmin_fdr_tmp),function(i){
    OG_pvals_tmp = OG_pvals[climate_var == picmin_fdr_tmp$climate_var[i] &
                              Orthogroup == picmin_fdr_tmp$Orthogroup[i],][order(epval_final),][1:picmin_fdr_tmp$config_est[i],]
    if(max(OG_pvals_tmp$epval_final > 0.1)){
      OG_pvals_tmp = OG_pvals_tmp[epval_final < 0.1,]
    }
    return(OG_pvals_tmp[,.(species,Orthogroup)])
  }))
  
  # Merge this with the GO terms to get species contributing to which GO terms...
  species_contributing_merge = merge(unique(species_contributing_tmp),OG_GO_tmp,by = "Orthogroup",allow.cartesian = T)
  
  if(length(unique(OG_GO_tmp$go_id)) > 2){
    
    # We now want to run semantic similarity clustering on the go terms to order them in a way that reflects semantic similarity
    semantic_dist = mgoSim(unique(species_contributing_merge$go_id), 
                           unique(species_contributing_merge$go_id), 
                           semData=AtGO, measure="Wang", combine=NULL)
    colnames(semantic_dist) <- rownames(semantic_dist) <- unique(species_contributing_merge$go_name)
    go_clusters = hclust(as.dist(1 - semantic_dist))
    # plot(go_clusters)
    
    # Plot these as a vertical dendrogram
    go_dendro <- as.dendrogram(go_clusters)
    # Rectangular lines
    go_dendro_plot <- dendro_data(go_dendro, type = "rectangle")
    go_dendro_fig <- ggplot(segment(go_dendro_plot)) + 
      geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
      coord_flip() + 
      theme_void() +
      scale_x_continuous(expand = c(0, 0),limits = c(0.5,max(go_dendro_plot$labels$x)+0.5)) 
    go_dendro_fig
    
    # Plot the GO-Species heatmap...
    go_order = data.table(go_dendro_plot$labels)
    go_order = go_order[order(x),]
    
  }
  
  # Make the final plottable output...
  species_GO_counts = species_contributing_merge[,.(`Orthogroup N` = nrow(.SD)),by = .(go_name,species)]
  if(length(unique(OG_GO_tmp$go_id)) > 2){
    species_GO_counts$go_name_F = factor(species_GO_counts$go_name,levels = go_order$label)
  } else {
    species_GO_counts$go_name_F = factor(species_GO_counts$go_name)
  }
  species_GO_counts$species_F = factor(species_GO_counts$species,levels = get_species_order()$species)
  species_GO_counts$climate_var = var_titles[x]
  species_GO_counts
})

# Plot these all together...
go_heat_data = rbindlist(GO_species_heat_list)
go_heat_data$climate_var_F = factor(go_heat_data$climate_var,levels = unique(go_heat_data$climate_var))
species_GO_heat = ggplot(go_heat_data,aes(y = go_name_F,x = species_F,fill = `Orthogroup N`)) +
  geom_tile() +
  scale_fill_viridis(option = 'D') +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "gray"),
        axis.title = element_blank(),
        title = element_text(size = 14),
        strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = 'gray20',colour = 'gray10')) +
  facet_grid(climate_var_F~.,scales = 'free_y', space = "free_y",switch = 'both') +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5))
species_GO_heat

# Add in bars for the significance of each of these...
for(i in 1:length(climate_enrichment)){
  climate_enrichment[[i]]$climate_var = var_titles[i]
}

go_bars_data = unique(go_heat_data[,.(go_name,go_name_F,climate_var,climate_var_F)])
for(i in 1:nrow(go_bars_data)){
  go_bars_data$fdr[i] = rbindlist(climate_enrichment)[climate_var == go_bars_data$climate_var_F[i] &
                                                        go_name == go_bars_data$go_name[i],fdr]
}
go_bars_fig = ggplot(go_bars_data,aes(y = go_name_F,x = -log10(fdr))) +
  geom_bar(stat = 'identity',colour = 'black',fill = 'gold2') +
  facet_grid(climate_var_F~.,scales = 'free_y', space = "free_y") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 14,angle = 45,hjust = 1),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_blank(),
        legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        strip.text = element_blank(),
        strip.background = element_blank()) +
  labs(x = expression(-log[10]~FDR)) +
  geom_vline(xintercept = -log10(0.1),linetype = 'dashed')
go_bars_fig

# Merge with the heatmap...
GO_heat_bars_merge = cowplot::plot_grid(species_GO_heat,
                                        go_bars_fig,
                                        ncol = 2,rel_widths = c(5,1),axis = 'tblr',align = 'h')


# Plot the precipitation network ------------------------------------------
clim_network = OG_string_codes[Orthogroup %in% climate_OG_convergent$precip_vars,STRING_id]
Athal_string_db$plot_network(clim_network)

# Fetch interactions to plot...
interactions_to_plot = data.table(Athal_string_db$get_interactions(clim_network))
interactions_to_plot = interactions_to_plot[combined_score > 400,]

library(GGally)
library(ggnet)
library(network)
library(sna)
library(ggplot2)

network_obj <- network(as.matrix(interactions_to_plot[,1:2]))
network_fig = ggnet2(network_obj)

# Loop over some choice GO terms and identify the nodes with genes related to that function...
go_vars = climate_enrichment$precip_vars$go_name
go_var_nodes = rbindlist(lapply(go_vars,function(go_var){
  OG_tmp = OG_GO[go_name == go_var & Orthogroup %in% climate_OG_convergent$precip_vars,Orthogroup]
  string_tmp = OG_string_codes[Orthogroup %in% OG_tmp,"STRING_id"]
  string_tmp$go_name = go_var
  string_tmp
}))

perGO_networks = lapply(unique(go_var_nodes$go_name),function(x){
  network_obj_tmp = network_obj
  network_obj_tmp %v% "GO1" = ifelse(sapply(network_obj$val,'[[',2) %in% go_var_nodes[go_name == x,STRING_id],"Yes",NA)

  
  set.seed(1001)
  ggnet2(network_obj_tmp, 
         color = "GO1",
         palette = c("Yes" = "orange2"),
         alpha = 0.75,
         size = 4) +
    theme(legend.position = "none") +
    ggtitle(gsub("Regulation Of ","Regulation Of\n",stringr::str_to_title(x))) 
})
names(perGO_networks) = unique(go_var_nodes$go_name)

network_egs = cowplot::plot_grid(perGO_networks$`response to auxin`,
                                 perGO_networks$`regulation of flower development`,
                                 perGO_networks$`regulation of stomatal movement`,
                                 perGO_networks$`root hair cell tip growth`,ncol = 2,align = 'hv',axis = 'tblr')
network_egs


# Plot Figure 3 -----------------------------------------------------------
pdf("figs/Figure3_network_interactions_and_GO.pdf",width = 14,height=8)
cowplot::plot_grid(
  cowplot::plot_grid(string_perm_res_dots,
                     network_egs,
                     ncol = 1,rel_heights = c(2,4),labels = c("A","B"),label_size = 32,label_y = c(1,1.05)),
  GO_heat_bars_merge,
  ncol = 2,rel_widths = c(1,1.5),labels = c("","C"),label_size = 32,align = "h",axis = "tblr")
dev.off()

# Also save the enrichment results...
enrichment_export = rbindlist(climate_enrichment)
write.table(enrichment_export,
            'tables/TableSX_fig3_GO_enrichment_results.tsv',
            row.names = F,quote = F,sep = '\t')
