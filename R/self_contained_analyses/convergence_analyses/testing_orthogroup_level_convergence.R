# Orthogroup-level association with environment...
lib <- c("ggtree","ape","VGAM","ggExtra","WeMix","pbmcapply","parallel","wCorr","data.table","ggplot2","viridis","ggridges","dplyr","patchwork","readr")
sapply(lib,library,character.only=T)

# Where are the blast results?
OG_dir <- "outputs/orthology/Results_210922_17_genomes_with_conifers_newPtrem_noAA_filter/Orthogroups/"
OG_files <- list.files(blast_dir,pattern="Blast")

# Where is all the metadata...
vcf_genome_map <- read.table("metadata/vcf_genome_gff_210922_map.txt",fill=T)
dataset_meta <- vcf_genome_map[,c("V8","V5")]
colnames(dataset_meta) <- c("dataset","genome")

# What are the OF codes for each of the genomes
genome_OF_codes <- data.frame(prot=c("Ahalleri.faa",
                                     "Alyrata.faa",
                                     "Athaliana.faa",
                                     "Atubercatus.faa",
                                     "Bstricta.faa",
                                     "Crubella.faa",
                                     "Egrandis.faa",
                                     "Hannuus.faa",
                                     "Mtruncatula.faa",
                                     "Pabies.faa",
                                     "Pdeltoides.faa",
                                     "Phallii.faa",
                                     "Pmenziesii.faa",
                                     "Ptaeda.faa",
                                     "Ptremula.faa",
                                     "Ptrichocarpa.faa",
                                     "Qpetraea.faa"))
genome_OF_codes$genome <- gsub(".faa","",genome_OF_codes$prot)
genome_OF_codes$OF_code <- (1:nrow(genome_OF_codes))-1

# Read in all of our liftovers from gff to OF codes
gff_OF_liftovers <- list.files("data/reference_genomes/",recursive = TRUE,pattern="_proteome_to_OF_id_map.txt")
gff_OF_liftovers <- lapply(genome_OF_codes$genome,function(x) grep(x,gff_OF_liftovers,value=T))
names(gff_OF_liftovers) <- genome_OF_codes$genome

# Fetch the sequence IDs...
OG_sequences <- read.table(paste0(dirname(OG_dir),"/WorkingDirectory/SequenceIDs.txt"),fill=T)[,1:2]
colnames(OG_sequences) <- c("OF_gene","gene_name")
OG_sequences$OF_gene <- gsub(":","",OG_sequences$OF_gene)
OG_sequences$OF_genome <- sapply(strsplit(OG_sequences$OF_gene,"_"),'[[',1)
for(genome in genome_OF_codes$genome){
  OF_code <- genome_OF_codes[genome_OF_codes$genome == genome,"OF_code"]
  OG_sequences[OG_sequences$OF_genome == as.character(OF_code),"gene_name"] <- paste0(genome,"_",OG_sequences[OG_sequences$OF_genome == as.character(OF_code),"gene_name"])
}

# Transform names so they match
OG_sequences$gene_name <- gsub("transcript:","transcript_",OG_sequences$gene_name)
OG_sequences$gene_name <- gsub("RNA:","RNA_",OG_sequences$gene_name)

# Summarise orthogroups ---------------------------------------------------
OG_gene_counts <- data.frame(fread(paste0(OG_dir,"Orthogroups.GeneCount.tsv")))

# How many orthogroups after removing any with >N in one species...
max_cutoffs <- c(2,5,10,20,30,40,50)
row_maxes <- apply(OG_gene_counts[,genome_OF_codes$genome],1,max)
row_empties <- apply(OG_gene_counts[,genome_OF_codes$genome],1,min)
zero_count <- function(x){ length(x[x==0]) }
row_zeroes <- apply(OG_gene_counts[,2:ncol(OG_gene_counts)],1,zero_count)

cutoff_props <- data.frame(cutoff = max_cutoffs,
                           prop_left = sapply(max_cutoffs,function(x) length(row_maxes[row_maxes <= x])/length(row_maxes)))

hist(OG_gene_counts$Total[OG_gene_counts$Total < 20] )
median(OG_gene_counts$Total)

ggplot(cutoff_props,aes(x=cutoff,y=prop_left))+
  geom_line()+
  geom_hline(yintercept = 1,linetype="dashed")+
  labs(y="Proportion of OG remaining",x="Max genes per genome")

# For now, let's take a subset of HQ orthogroups...
OG_HQ <- which(row_zeroes <= 7 & row_maxes < 10)
length(OG_HQ)/nrow(OG_gene_counts)
OG_HQ_subset <- OG_gene_counts$Orthogroup[OG_HQ]

# Prepare GEA data -------------------------------------------------------
focal_datasets <- list.files("outputs/GEA_res")
focal_climate <- "mean_temp"

# Remove any focal_datasets without full wza outputs...
focal_datasets <- na.omit(sapply(focal_datasets,function(x){
  if(file.exists(paste0("outputs/GEA_res/",x,"/",focal_climate,"_WZA_TC_allgenes.rds"))){
    return(x)
  } else {
    return(NA)
  }
}))

# Read in all of our GEA results
all_gea_res <- pbmclapply(focal_datasets,function(dataset){
  
  # Get raw GEA
  wza <- readRDS(paste0("outputs/GEA_res/",dataset,"/",focal_climate,"_WZA_TC_allgenes.rds"))
  wza$gea_gene <- wza$gene_id
  
  # Merge these with the relevant OF_code...
  genome <- na.omit(dataset_meta[dataset_meta$dataset == dataset,"genome"])
  OF_codes <- read.table(paste0("data/reference_genomes/",gff_OF_liftovers[genome]),header=T)
  if(colnames(OF_codes)[3]=="seqid"){
    OF_codes$gea_gene <- paste0(OF_codes$seqid,":",OF_codes$start,"-",OF_codes$end)
  }
  OF_codes$full_OF <- paste0(genome_OF_codes[genome_OF_codes$genome == genome,"OF_code"],"_",OF_codes$OF_ID)
  wza_merge <- merge(wza[,c("gea_gene","weiZ")],OF_codes[,c("gea_gene","full_OF")])
  
  # Mark
  wza_merge$dataset <- dataset
  wza_merge$genome <- genome
  return(wza_merge)
},mc.cores=6)
names(all_gea_res) <- focal_datasets


# OG-level GEA ------------------------------------------------------------
test_OG <- OG_HQ_subset[1]

# Fetch the tree
OG_tree <- read.tree(paste0(dirname(OG_dir),"/Resolved_Gene_Trees/Resolved_Gene_Trees/",test_OG,"_tree.txt"))

# Build tree metadata
OG_tree_metadata <- data.frame(gene_name=OG_tree$tip.label,
                               tip.pos=1:length(OG_tree$tip.label),
                               OG=test_OG)
OG_tree_metadata <- merge(OG_tree_metadata,OG_sequences[,c("OF_gene","gene_name")],by="gene_name")

# Fetch all of the relevant GEA results...
OG_GEA_res <- data.frame(rbindlist(lapply(all_gea_res,function(gea){
  
  # Subset for OG genes
  gea_sub <- gea[gea$full_OF %in% OG_tree_metadata$OF_gene,]
  
  if(nrow(gea_sub)==0){
    return(NULL)
  } else {
    gea_sub$OF_gene <- gea_sub$full_OF
    gea_sub <- merge(gea_sub,OG_tree_metadata[,c("tip.pos","OF_gene")],by="OF_gene")
    return(gea_sub)
  }
  
})))

# Identify which tips are duplicated and split these...
duplicated_tips <- table(OG_GEA_res$OF_gene)

# Drop any tips that we don't have data for...
to_drop <- OG_tree_metadata$OF_gene[!(OG_tree_metadata$OF_gene %in% OG_GEA_res$OF_gene)]
dropped_OG_tree <- drop.tip(OG_tree,unique(OG_tree_metadata[OG_tree_metadata$OF_gene %in% to_drop,"tip.pos"]))

# Build tree metadata
OG_tree_metadata <- data.frame(gene_name=dropped_OG_tree$tip.label,
                               tip.pos=1:length(dropped_OG_tree$tip.label),
                               OG=test_OG)
OG_tree_metadata <- merge(OG_tree_metadata,OG_sequences[,c("OF_gene","gene_name")],by="gene_name")

# Build a GEA matrix to plot alongside the
gea_matrix <- matrix(nrow=nrow(OG_tree_metadata),ncol=max(duplicated_tips))
rownames(gea_matrix) <- dropped_OG_tree$tip.label
for(i in 1:nrow(gea_matrix)){
  OF_gene_tmp <- OG_tree_metadata[OG_tree_metadata$gene_name==rownames(gea_matrix)[i],"OF_gene"]
  for(j in 1:duplicated_tips[OF_gene_tmp]){
    gea_matrix[i,j] <- OG_GEA_res[OG_GEA_res$OF_gene == OF_gene_tmp,"weiZ"][j]
  }
}
# gea_matrix[is.na(gea_matrix)] <- 0
colnames(gea_matrix) <- paste0("V",1:ncol(gea_matrix))

# Make a tree with phenos...
base_tree <- ggtree(dropped_OG_tree,ladderize = F,branch.length = "none")+ 
  theme_tree()+
  xlim(0, 2)


# Add to tree
final_wza_tree <- gheatmap(base_tree, gea_matrix,colnames=FALSE, legend_title="GEA WZA")+
  scale_fill_viridis(option="A")+
  scale_x_ggtree()+
  theme(legend.position = "top",
        legend.title = element_text(size=14),
        legend.text = element_text(size=13))+
  labs(fill="GEA WZA")

# # Build a distance matrix between orthologues...
# dropped_OG_tree_BL <- compute.brlen(dropped_OG_tree)

# Plot histograms...
# All wza
all_wza_histogram <- ggplot(OG_GEA_res,aes(x=weiZ))+
 # geom_histogram()+
  stat_density(fill="skyblue",alpha=0.6)+
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = median(OG_GEA_res$weiZ),colour="red2")+
  theme_minimal()+
  labs(y="Count",x="WZA")+
  ggtitle(paste0(test_OG," all WZA"))

# Best wza
OG_GEA_res_max <- data.frame(OG_GEA_res %>% group_by(dataset) %>% summarise(max_wza=max(weiZ)))
max_wza_histogram <- ggplot(OG_GEA_res_max,aes(x=max_wza))+
  #geom_histogram()+
  stat_density(fill="skyblue",alpha=0.6)+
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = median(OG_GEA_res_max$max_wza),colour="red2")+
  theme_minimal()+
  labs(y="Count",x="WZA (max)")+
  ggtitle(paste0(test_OG," max WZA per dataset"))

# Add all together
library(patchwork)
final_wza_tree|(all_wza_histogram/max_wza_histogram)





