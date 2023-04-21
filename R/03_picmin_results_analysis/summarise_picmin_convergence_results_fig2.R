# Analyse and summarise PicMin outputs...
lib <- c("biomaRt","org.At.tair.db","ggridges","ggrepel","cowplot","ggplot2","data.table","tidyr","pbmcapply","ggtree","ape","viridis","STRINGdb","qvalue","dplyr")
sapply(lib,library,character.only=T)
n_cores = 6
source("R/repadapt_functions.R")

# Which run name are we looking at?
run_name = "230321"
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"

# Set up the Phylo-Ordered datasets... ------------------------------------
# Fetch our Orthofinder phylogeny...
OF_tree <- read.tree("outputs/orthology/Results_221213_18_genomes_Ptaeda_isoforms_removed/Species_Tree/SpeciesTree_rooted.txt")

# Fetch our genome dataset map
dataset_meta <- read.table("metadata/vcf_genome_gff_220830_map.txt")
dataset_meta <- na.omit(dataset_meta[,c(8,5)])
colnames(dataset_meta) <- c("dataset","genome")

# Reorder metadata...
dataset_meta_ordered <- NULL
for(genome in OF_tree$tip.label){
  dataset_meta_ordered <- rbind(dataset_meta_ordered,dataset_meta[dataset_meta$genome == genome,])
}

# Reduce down for species...
dataset_meta_ordered$species = paste0(sapply(strsplit(dataset_meta_ordered$dataset,'_'),'[[',1),
                                      " ",
                                      sapply(strsplit(dataset_meta_ordered$dataset,'_'),'[[',2))


species_meta_ordered = unique(dataset_meta_ordered[,c("genome","species")])

# Fetch all of the climate clines -----------------------------------------
dataset_dirs <- list.files("outputs/GEA_res/",pattern = run_name)
dataset_dirs <- grep(".rds",invert=T,dataset_dirs,value = T)

climate_clines <- lapply(paste0("outputs/GEA_res/",dataset_dirs,"/climate_cline.tsv"),read.table,header=T)
names(climate_clines) <- gsub(paste0(run_name,"_"),"",dataset_dirs)

# Prepare Athal vs OG identifiers... --------------------------------------
if(!file.exists(paste0("data/OG_map_Athal_",output_name,".rds"))){
  OG_map <- readRDS(paste0("data/OF_OG_gea_gene_map_",output_name,".rds"))
  # OG_map <- data.table(readRDS("metadata/OF_to_genome_map.rds"))
  
  # Fetch all the headers in the proteome and convert to a data.table
  prot_path = "data/proteomes/Athaliana.faa"
  prot_headers <- system(paste0("grep '>' ",prot_path),intern = T)
  
  # Make into something new...
  Athal_prot_dd <- data.table(gene_name_noGenome = sapply(strsplit(prot_headers," gene="),'[[',1),
                              TAIR_gene = sapply(strsplit(prot_headers,"gene=| name="),'[[',2),
                              true_gene = sapply(strsplit(prot_headers,"name=| seq_id"),'[[',2),
                              OF_ID = 0:(length(prot_headers) - 1))
  Athal_prot_dd$gene_name_noGenome <- gsub(">","",Athal_prot_dd$gene_name_noGenome)
  Athal_prot_dd$TAIR_gene <- gsub("gene-","",Athal_prot_dd$TAIR_gene)
  
  # Now merge with teh OG_map
  OG_map_Athal <- OG_map[OG_map$genome == "Athaliana",]
  OG_map_Athal <- merge(OG_map_Athal,Athal_prot_dd,by="OF_ID")
  
  saveRDS(OG_map_Athal,paste0("data/OG_map_Athal_",output_name,".rds"))
} else {
  OG_map_Athal = data.table(readRDS(paste0("data/OG_map_Athal_",output_name,".rds")))
}

# And repeat for Mtrunc
if(!file.exists(paste0("data/OG_map_Mtrunc_",output_name,".rds"))){
  OG_map <- readRDS(paste0("data/OF_OG_gea_gene_map_",output_name,".rds"))
  # Fetch all the headers in the proteome and convert to a data.table
  prot_path = "data/proteomes/Mtruncatula.faa"
  prot_headers <- system(paste0("grep '>' ",prot_path),intern = T)
  prot_headers <- strsplit(prot_headers," ")
  
  # Make into something new...
  Mtrunc_prot_dd <- data.table(gene_name_noGenome = sapply(prot_headers,'[[',1),
                               gene_name = sapply(prot_headers,'[[',2),
                               OF_ID = 0:(length(prot_headers) - 1))
  Mtrunc_prot_dd$gene_name_noGenome <- gsub(">","",Mtrunc_prot_dd$gene_name_noGenome)
  Mtrunc_prot_dd$gene_name <- gsub("gene=","",Mtrunc_prot_dd$gene_name)
  
  # Now merge with teh OG_map
  OG_map_Mtrunc <- OG_map[OG_map$genome == "Mtruncatula",][order(as.integer(OF_ID)),]
  OG_map_Mtrunc <- merge(OG_map_Mtrunc,Mtrunc_prot_dd,by="OF_ID")
  
  # Convert to biomart names...
  OG_map_Mtrunc$biomart_gene = gsub("Medtr","MTR_",OG_map_Mtrunc$gene_name)
  
  saveRDS(OG_map_Mtrunc,paste0("data/OG_map_Mtrunc_",output_name,".rds"))
} else {
  OG_map_Mtrunc = data.table(readRDS(paste0("data/OG_map_Mtrunc_",output_name,".rds")))
}

# Fetch data and analyse --------------------------------------------------
picmin_outputs = readRDS(paste0("outputs/",output_name,"_picmin_results_doubleDS_pvals.rds"))
picmin_fdr = rbindlist(lapply(picmin_outputs,'[[',1))
OG_pvals = rbindlist(lapply(picmin_outputs,'[[',2))
# OG_pvals = OG_pvals[Orthogroup %in% picmin_fdr$Orthogroup,]

# Number of genes covered per species...
ngenes_per_species = OG_pvals[Orthogroup %in% unique(picmin_fdr$Orthogroup) & climate_var == "annual_precip",
                              .(Ngenes = sum(Ngenes_per_species)),by = species]
species_counts = table(dataset_meta_ordered$species)
for(species_tmp in names(species_counts[species_counts > 1])){
  ngenes_per_species[species == species_tmp,Ngenes := Ngenes/(species_counts[species])]
}

# Significance bars -------------------------------------------------------
# Plot stacked bars for convergence N...
picmin_fdr$fdr_level = ">0.5"
for(fdr in seq(0.5,0.1,-0.1)){
  picmin_fdr[picmin_fdr < fdr,fdr_level := paste0("<",fdr)]
}
picmin_fdr$clim_lab = stringr::str_to_title(gsub("_"," ",picmin_fdr$climate_var))
picmin_fdr$clim_lab_F = factor(picmin_fdr$clim_lab,levels = names(sort(table(picmin_fdr[picmin_fdr < 0.5,clim_lab]))))
picmin_fdr$fdr_level_F = factor(picmin_fdr$fdr_level,levels = paste0("<",seq(0.5,0.1,-0.1)))

fdr_bars = picmin_fdr[picmin_fdr < 0.5,.(OGn = nrow(.SD)),by = .(clim_lab_F,fdr_level_F)] |>
  ggplot(aes(y=clim_lab_F,x=OGn,fill=fdr_level_F))+
  geom_bar(stat="identity",show.legend = T)+
  theme_minimal()+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=18),
        strip.text = element_text(size=20),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        # legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        # legend.margin = margin(6, 6, 6, 6),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.key = element_rect(fill = "white", colour = "black"),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 14))+
  labs(x="Number of RAOs", fill = "FDR")+
  scale_fill_viridis_d(option = "C")

# Get some stats here
picmin_fdr[picmin_fdr < 0.5,.(total_fdr = nrow(.SD)),by = clim_lab_F][order(-total_fdr),]
picmin_fdr[picmin_fdr < 0.1,.(total_fdr = nrow(.SD)),by = clim_lab_F][order(-total_fdr),]
length(unique(picmin_fdr[picmin_fdr < 0.5,Orthogroup]))
length(unique(picmin_fdr[picmin_fdr < 0.2,Orthogroup]))

# Which genes appear most often?
picmin_fdr_signif <- data.frame(picmin_fdr)[picmin_fdr$picmin_fdr < 0.5,]
signif_OG_counts <- data.table(table(picmin_fdr_signif$Orthogroup))[order(-N),]

# Order these by most common...
multiple_picmin_OG = rbindlist(lapply(signif_OG_counts[N > 1,V1],function(OG){
  picmin_fdr_signif[picmin_fdr_signif$Orthogroup == OG,]
}))
multiple_picmin_OG


# Figure of Orthogroup with most hits -------------------------------------
OG_to_plot = signif_OG_counts$V1[1]

# Fetch all of the pvals for all climate vars for this orthogroups...
OG_pvals_top = OG_pvals[Orthogroup == OG_to_plot,]
OG_pvals_top$climate_labs = stringr::str_to_title(gsub("_"," ",OG_pvals_top$climate))

# Order these by picmin significance...
top_OG_picmin_res = picmin_fdr[Orthogroup == OG_to_plot,][order(picmin_fdr),.(picmin_p_adj,picmin_fdr,climate_var)]
top_OG_picmin_res$climate_lab = stringr::str_to_title(gsub("_"," ",top_OG_picmin_res$climate))
OG_pvals_top$climate_labs_F = factor(OG_pvals_top$climate_labs,levels = top_OG_picmin_res$climate_lab)

# Plot heatmap of species pvals for OG
OG_pvals_top$species_F = factor(OG_pvals_top$species,levels = get_species_order()[,2])
top_OG_heatmap = ggplot(OG_pvals_top,aes(y = climate_labs_F,x = species_F,fill = -log10(epval_final))) +
  geom_tile() +
  scale_fill_viridis(option = 'A') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),
        legend.position = "bottom",
        axis.title = element_blank(),
        title = element_text(size = 18)) +
  labs(fill = 'GEA\nmin epvalue') +
  scale_y_discrete(limits=rev) +
  ggtitle(OG_to_plot)

# Plot bars alongside
top_OG_picmin_res$climate_lab_F = factor(top_OG_picmin_res$climate_lab,levels = top_OG_picmin_res$climate_lab)
top_OG_bars = ggplot(top_OG_picmin_res,aes(y = climate_lab_F,x = -log10(picmin_fdr))) +
  geom_bar(stat = 'identity') +
  theme_minimal() +
  scale_y_discrete(limits=rev) +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        panel.grid.major.y = element_blank()) +
  labs(x = '-log10(PicMin FDR)') +
  geom_vline(xintercept = -log10(0.5),colour = "red2")

# Plot all of the individual pvalue distributions, in order of significance...
top_OG_pval_distributions = ggplot(OG_pvals_top,aes(x = epval_final)) +
  geom_histogram(bins = 20) +
  facet_wrap(~climate_labs_F,nrow = 3) +
  theme_bw() +
  theme(strip.text = element_text(size = 10),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14)) +
  labs(x = "GEA Minimum P-value (Paralog-corrected)",y = "Count")

# Combine
top_OG_combine = cowplot::plot_grid(
  cowplot::plot_grid(top_OG_heatmap,top_OG_bars,
                     ncol = 2,rel_widths = c(5,1),
                     align = "h",axis = "tblr"),
  top_OG_pval_distributions,labels = "AUTO",label_size = 32,ncol = 1,align = "v",axis = "tblr",rel_heights = c(2,1),vjust = c(1.5,0))
pdf("figs/FigureSX_picmin_topOG_summary_heatmap.pdf",width=13,height=10)
top_OG_combine
dev.off()

# Figure 2 heatmaps -------------------------------------------------------
picmin_fdr_signif = picmin_fdr[picmin_fdr < 0.5,]
# Fetch all the datasets associated with each signif OG i.e. a pvalue < 0.1
picmin_fdr_signif_dataset <- rbindlist(lapply(1:nrow(picmin_fdr_signif),function(x){
  
  # Find the original pvals and take 1:config.est
  tmp_gea <-  OG_pvals[climate_var == picmin_fdr_signif$climate_var[x] &
                         Orthogroup == picmin_fdr_signif$Orthogroup[x],][order(epval_final),][1:picmin_fdr_signif$config_est[x]]
  # If the cutoff pvalues is >0.1, only take those with pval <0.1
  if(max(tmp_gea$epval_final) > 0.1){
    tmp_gea = tmp_gea[epval_final < 0.1,]
  }
  tmp_gea$climate_var <- picmin_fdr_signif$climate_var[x]
  
  if("species" %in% colnames(tmp_gea)){
    tmp_gea$dataset = tmp_gea$species
  }
  return(tmp_gea[,c("Orthogroup","climate_var","dataset")])
}))


# Dataset vs Climate variable...
picmin_fdr_signif_dataset$dataset_clim <- paste0(picmin_fdr_signif_dataset$dataset,":",picmin_fdr_signif_dataset$climate_var)
dataset_clim_counts <- data.frame(table(picmin_fdr_signif_dataset$dataset_clim)) %>%
  separate(col = "Var1",into=c("species","climate"),sep=":")

# We also want to get these as proportions of those tested among each species/climate...
dataset_clim_tested = OG_pvals[,.(testedN = nrow(.SD)),by = .(species,climate_var)]

# Set factors for plotting...
dataset_clim_counts$species_F <- factor(dataset_clim_counts$species,levels=species_meta_ordered$species)
dataset_clim_counts$clim_F <- factor(stringr::str_to_title(gsub("_"," ",dataset_clim_counts$climate)),levels=levels(picmin_fdr$clim_lab_F))

# Merge these...
dataset_clim_counts = merge(dataset_clim_counts,dataset_clim_tested)

# Plot climate vs species heatmap
climate_species_heat = ggplot(dataset_clim_counts,aes(y=species_F,x=clim_F,fill=Freq/testedN))+
  geom_tile()+
  scale_fill_viridis(option="A")+
  theme_minimal()+
  theme(axis.title=element_blank(),
        axis.text.x = element_text(size=10,angle=45,hjust=1),
        panel.grid=element_blank(),
        panel.background = element_rect(fill = "gray"))+
  labs(fill = expression(N/N[Tested])) +
  ggtitle("Contribution to RAOs (FDR < 0.5)\nby dataset and climate")
climate_species_heat

data.table(dataset_clim_counts)[,.(all_contributions = mean(.SD$Freq/.SD$testedN)),by = species][order(-all_contributions),]


# Dataset vs Genes variable...
picmin_fdr_signif_dataset$dataset_OG <- paste0(picmin_fdr_signif_dataset$dataset,":",picmin_fdr_signif_dataset$Orthogroup)
dataset_OG_counts <- data.table(table(picmin_fdr_signif_dataset$dataset_OG)) %>%
  separate(col = "V1",into=c("dataset","Orthogroup"),sep=":")
dataset_OG_counts = data.table(dataset_OG_counts)
dataset_OG_counts$dataset <- gsub(paste0(run_name,"_"),"",dataset_OG_counts$dataset)

# Limit to OG where fdr < 0.5 in more than one variable...
OG_to_plot = picmin_fdr[picmin_fdr < 0.5,Orthogroup]
# OG_to_plot = unique(OG_to_plot[duplicated(OG_to_plot)])

# Organise species by phylogeny and organise OG by how common they are...
dataset_OG_counts$dataset_F = factor(dataset_OG_counts$dataset,levels = unique(dataset_meta_ordered$species))
dataset_OG_counts$Orthogroup_F = factor(dataset_OG_counts$Orthogroup,levels = rev(data.table(table(dataset_OG_counts$Orthogroup))[order(-N),V1]))

# Now replace Orthogroups with gene names
dataset_OG_counts$TAIR_genes = sapply(dataset_OG_counts$Orthogroup,function(OG){
  paste(OG_map_Athal[Orthogroup == OG,true_gene],collapse = ' / ')
})
dataset_OG_counts$TAIR_genes_F = factor(dataset_OG_counts$TAIR_genes,levels = rev(data.table(table(dataset_OG_counts$TAIR_genes))[order(-N),V1]))

# Filter again
OG_species_counts = table(unique(dataset_OG_counts[,.(dataset,Orthogroup)])$Orthogroup)
# Let's only plot genes that have contributions from at least 5 species...
OG_to_plot = names(OG_species_counts[OG_species_counts >= 5])
OG_species_heat = ggplot(dataset_OG_counts[Orthogroup %in% OG_to_plot,],aes(x=dataset_F,y=TAIR_genes_F,fill=N))+
  geom_tile()+
  scale_fill_viridis(option="A") +
  theme_minimal() +
  theme(axis.title=element_blank(),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        panel.background = element_rect(fill = "gray"),
        panel.grid = element_blank())+
  labs(fill = "N Climate\nVars") +
  geom_vline(xintercept = 5.5,size = 2, colour = 'green3',linetype = 'solid')
OG_species_heat

# Also produce a bar graph showing the number of times an Orthgroup pops up and the total number of orthogroups
OG_count_bar_data = dataset_OG_counts[,.(`Species N` = nrow(.SD),
                                         `Total N` = sum(.SD$N) - nrow(.SD)), by = Orthogroup]
OG_count_bar_data = merge(OG_count_bar_data,unique(dataset_OG_counts[,.(Orthogroup,TAIR_genes)]),by = "Orthogroup")
OG_count_bar_data$TAIR_genes_F = factor(OG_count_bar_data$TAIR_genes,levels = levels(dataset_OG_counts$TAIR_genes_F))
OG_species_bars = melt(OG_count_bar_data[Orthogroup %in% OG_to_plot,]) |>
  mutate(variable_F = factor(variable,levels = c("Total N","Species N"))) |>
  ggplot(aes(y=TAIR_genes_F,x=value,fill=variable_F))+
  geom_bar(stat="identity",show.legend = T)+
  theme_minimal()+
  theme(axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=18),
        strip.text = element_text(size=20),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        # legend.box.just = "right",
        legend.box.background = element_rect(color="black", size=1),
        # legend.margin = margin(6, 6, 6, 6),
        legend.position = c(.95, .05),
        legend.justification = c("right", "bottom"),
        legend.key = element_rect(fill = "white", colour = "black"),
        # legend.title = element_text(size = 16),
        legend.text = element_text(size = 14),
        legend.title = element_blank())+
  labs(x="Contributions to repeatability\n(FDR < 0.5)")+
  scale_fill_brewer(palette = "Set1")

OG_species_combined = cowplot::plot_grid(OG_species_heat + 
                                           ggtitle("Contribution of species to most\nsignificant (FDR < 0.5) RAOs") + 
                                           theme(legend.position = "bottom",title = element_text(size = 14)),
                                         OG_species_bars + 
                                           theme(axis.text.y = element_blank()) + 
                                           xlab("Sum"),
                                         ncol = 2,align="h",axis = "tblr",rel_widths = c(2,0.75))
OG_species_combined

# We also want to know how many times both species were tested together...
# Make a list of all climate-OG that were tested per species
species_tested = pbmclapply(unique(OG_pvals$species),function(spp) apply(OG_pvals[species == spp,.(climate_var,Orthogroup)],1,paste,collapse = ":"),
                            mc.cores = n_cores)
names(species_tested) = unique(OG_pvals$species)
species_tested_overlap = crossprod(table(stack(species_tested)))
diag(species_tested_overlap) = NA
pheatmap::pheatmap(species_tested_overlap)
species_tested_overlap_melt = melt(species_tested_overlap)
colnames(species_tested_overlap_melt) = c("Var1","Var2","testedN")

# Species vs Species - Same Orthogroup, same climate variable
picmin_fdr_signif_dataset$clim_OG <- paste0(picmin_fdr_signif_dataset$climate,":",picmin_fdr_signif_dataset$Orthogroup)
species_species_mat <- matrix(nrow=length(unique(picmin_fdr_signif_dataset$dataset)),
                              ncol=length(unique(picmin_fdr_signif_dataset$dataset)))
for(i in 1:ncol(species_species_mat)){
  for(j in 1:nrow(species_species_mat)){
    if(i > j){
      tmp <- picmin_fdr_signif_dataset[picmin_fdr_signif_dataset$dataset %in% unique(picmin_fdr_signif_dataset$dataset)[c(i,j)]]
      clim_OG_counts <- table(tmp$clim_OG)
      observed_overlap =  length(clim_OG_counts[clim_OG_counts > 1])
      tmp_testedN = species_tested_overlap[unique(picmin_fdr_signif_dataset$dataset)[i],
                                           unique(picmin_fdr_signif_dataset$dataset)[j]]
      
      
      species_species_mat[i,j] <- species_species_mat[j,i] <- observed_overlap / tmp_testedN
    }
  }
}
colnames(species_species_mat) <- rownames(species_species_mat) <- unique(picmin_fdr_signif_dataset$dataset)

# Climate variable vs Climate variable
# How many Orthogroups are called that match across climate vars
perClimate_convergent_OG = lapply(unique(picmin_fdr$climate_var),function(clim) picmin_fdr[picmin_fdr < 0.5 & climate_var == clim,Orthogroup])
names(perClimate_convergent_OG) = unique(picmin_fdr$climate_var)
# Intersect list
climate_OG_overlap = crossprod(table(stack(perClimate_convergent_OG)))
diag(climate_OG_overlap) = NA
pheatmap::pheatmap(climate_OG_overlap)

# Save as a supp fig
pdf("figs/FigureSX_climate_climate_overlap_heatmap.pdf",width = 8,height = 8)
print(pheatmap::pheatmap(climate_OG_overlap))
dev.off()

to_plot <- melt(species_species_mat)
to_plot$Var1_F <- factor(to_plot$Var1,levels=species_meta_ordered$species)
to_plot$Var2_F <- factor(to_plot$Var2,levels=species_meta_ordered$species)


# Plot the dataset vs dataset matrix
to_plot_merge = data.table(merge(to_plot,species_tested_overlap_melt,all.x = T))
to_plot_merge[Var1 == Var2, testedN := 1]
dataset_dataset_heatmap = ggplot(to_plot_merge,aes(x=Var2_F,y=Var1_F,fill=value/testedN))+
  geom_tile()+
  scale_fill_viridis(option="A")+
  theme_minimal()+
  theme(axis.title=element_blank(),
        axis.text.x = element_text(size=8,angle=45,hjust=1),
        panel.grid = element_blank())+
  ggtitle("Shared contribution to RAOs (FDR < 0.5)\n(same orthogroup, same climate)")+
  labs(fill=expression(N/N[Tested]))
dataset_dataset_heatmap
species_tested_overlap_melt

# Combine all convergence result plots... ---------------------------------
pdf("figs/Figure2_convergence_OG_results.pdf",width = 22,height = 11)
cowplot::plot_grid(fdr_bars,
                   cowplot::plot_grid(climate_species_heat,dataset_dataset_heatmap,ncol=1,align = "v",axis = "tblr",labels = c("B","C"),label_size = 32),
                   OG_species_combined,
                   ncol = 3,labels = c("A","","D"),label_size = 32,label_x = c(0,0,0.1),rel_widths = c(0.6,0.8,1.2))
dev.off()


# Build a summary table of strongest Orthogroups ---------------
fdr_cutoff = 0.5
most_significant_OG = unique(picmin_fdr[picmin_fdr < fdr_cutoff,Orthogroup])
most_significant_OG = OG_count_bar_data[Orthogroup %in% most_significant_OG,][order(-`Species N`),Orthogroup]

signif_OG_summary <- rbindlist(lapply(most_significant_OG,function(OG){
  clim_tmp = paste(stringr::str_to_title(gsub("_"," ",picmin_fdr[Orthogroup == OG & picmin_fdr < 0.5,climate_var])),collapse = "/")
  true_gene_tmp = paste(OG_map_Athal[Orthogroup == OG,true_gene],collapse = "/")
  OG_count = rowSums(OG_count_bar_data[Orthogroup == OG,.(`Species N`,`Total N`)])
  dataset_count_tmp = OG_count_bar_data[Orthogroup == OG,`Species N`]
  out = data.table(Orthogroup = OG,
                   convergence_N = OG_count,
                   dataset_N = dataset_count_tmp,
                   climate_vars = clim_tmp,
                   Genes = true_gene_tmp)
}))

# Prettify
signif_OG_summary = signif_OG_summary[order(-dataset_N,-convergence_N),]
colnames(signif_OG_summary) = c("Orthogroup","Total Convergent N (FDR < 0.5)","Species Convergent N (FDR < 0.5)","Climate Variables","Athaliana Genes")

# Save
write.csv(signif_OG_summary,"tables/TableSX_dataset_OG_heatmap_annotated_genes.csv",row.names = F)


# Another version of this with the TAIR descriptions directly -------------
# And make table...
fdr05_OG_summary <- rbindlist(lapply(signif_OG_summary$Orthogroup,function(OG){
  clim_tmp = paste(stringr::str_to_title(gsub("_"," ",picmin_fdr[Orthogroup == OG & picmin_fdr < 0.5,climate_var])),collapse = "/")
  true_gene_tmp = OG_map_Athal[Orthogroup == OG,true_gene]
  TAIR_gene_tmp = OG_map_Athal[Orthogroup == OG,TAIR_gene]
  
  out = data.table(Orthogroup = OG,
                   `Climate Variable` = clim_tmp,
                   `Athaliana Genes` = true_gene_tmp,
                   `TAIR Genes` = TAIR_gene_tmp)
  if(nrow(out) > 1){
    out[2:nrow(out),1:2] = ""
  }
  out
}))

write.csv(fdr05_OG_summary,"tables/TableSX_fdr05_TAIR_annotated_genes.csv",row.names = F)



# Add to these biomart function -------------------------------------------
# Add in GO annotations
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

# Do the stupid gene seps...
go_res = rbindlist(lapply(signif_OG_summary$Orthogroup,function(OG){
  mart_tmp = athal_universe_OG[Orthogroup == OG,]
  out = data.table(unique(mart_tmp[,.(go_id,go_name)]),
                   Orthogroup = OG)
}))

# Group all of these together into another supp table...
most_common_go = rev(sort(table(go_res$go_name)))
signif_OG_summary_functions = rbindlist(lapply(names(most_common_go),function(go){
  tmp = go_res[go_name == go,]
  merge(tmp,signif_OG_summary[,.(Orthogroup,`Athaliana Genes`)])
}))

# Save
write.csv(signif_OG_summary_functions[go_id != "",],"tables/TableSX_dataset_OG_heatmap_annotated_genes_GO_functions.csv",row.names = F)

# Build a summary of climate change genes ---------------------------------
fdr_cutoff = 0.5
climchange_OG = unique(picmin_fdr[picmin_fdr < fdr_cutoff &
                                    climate_var %in% c("tmax_clim_change","prec_clim_change"),.(Orthogroup,climate_var,picmin_fdr)])
climchange_OG = climchange_OG[order(picmin_fdr),]

# And make table...
climchange_OG_summary <- rbindlist(lapply(climchange_OG$Orthogroup,function(OG){
  clim_tmp = paste(stringr::str_to_title(gsub("_"," ",picmin_fdr[Orthogroup == OG & 
                                                                   picmin_fdr < 0.5,climate_var])),collapse = "/")
  true_gene_tmp = OG_map_Athal[Orthogroup == OG,true_gene]
  TAIR_gene_tmp = OG_map_Athal[Orthogroup == OG,TAIR_gene]
  fdr_tmp = climchange_OG[Orthogroup == OG,picmin_fdr]
  
  
  # # And get the GO terms and stitch
  # go_tmp = paste(athal_universe_OG[Orthogroup == OG & go_name != "",go_name],collapse = '/')
  
  out = data.table(Orthogroup = OG,
                   `Climate Variable` = clim_tmp,
                   `PicMin FDR` = fdr_tmp,
                   `Athaliana Genes` = true_gene_tmp,
                   `TAIR Genes` = TAIR_gene_tmp)
  if(nrow(out) > 1){
    out[2:nrow(out),1:3] = ""
  }
  out
}))

# Prettify
write.csv(climchange_OG_summary,"tables/TableSX_climchange_annotated_genes.csv",row.names = F)
