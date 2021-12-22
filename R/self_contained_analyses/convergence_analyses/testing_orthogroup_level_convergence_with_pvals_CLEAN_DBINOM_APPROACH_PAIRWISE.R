# Orthogroup-level association with environment...
lib <- c("mvmeta","qvalue","tidyr","ggtree","ape","VGAM","ggExtra","WeMix","pbmcapply","parallel","wCorr","data.table","ggplot2","viridis","ggridges","dplyr","patchwork","readr")
sapply(lib,library,character.only=T)

# Function Library --------------------------------------------------------
# Fetch for pairwise_c_hyper
# devtools::install_github("samyeaman/dgconstraint", build_vignettes = TRUE)
library(dgconstraint)

############################################################################
focal_climate <- "mean_temp_cold_quarter"

n_cores=6
orthogroup_cutoff <- 20

# Where are the blast results?
OG_dir <- "outputs/orthology/Results_210922_17_genomes_with_conifers_newPtrem_noAA_filter/Orthogroups/"

# Where is all the metadata...
vcf_genome_map <- read.table("metadata/vcf_genome_gff_210922_map.txt",fill=T)
dataset_meta <- vcf_genome_map[,c("V8","V5")]
colnames(dataset_meta) <- c("dataset","genome")

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

# Remove empties
OG_sequences <- OG_sequences[OG_sequences$V2 != "",]

# Format text
colnames(OG_sequences) <- c("OF_gene","gene_name")
OG_sequences$OF_gene <- gsub(":","",OG_sequences$OF_gene)
OG_sequences$OF_genome <- sapply(strsplit(OG_sequences$OF_gene,"_"),'[[',1)
OG_sequences$gene_name_noGenome <- OG_sequences$gene_name
for(genome in genome_OF_codes$genome){
  OF_code <- genome_OF_codes[genome_OF_codes$genome == genome,"OF_code"]
  OG_sequences[OG_sequences$OF_genome == as.character(OF_code),"gene_name"] <- paste0(genome,"_",OG_sequences[OG_sequences$OF_genome == as.character(OF_code),"gene_name"])
}

# Transform names so they match
OG_sequences$gene_name <- gsub("transcript:","transcript_",OG_sequences$gene_name)
OG_sequences$gene_name_noGenome <- gsub("transcript:","transcript_",OG_sequences$gene_name_noGenome)

OG_sequences$gene_name <- gsub("RNA:","RNA_",OG_sequences$gene_name)
OG_sequences$gene_name_noGenome <- gsub("RNA:","RNA_",OG_sequences$gene_name_noGenome)

# We also now want to attach the Orthogroup to each OG_sequence...
all_OG <- data.frame(fread(paste0(OG_dir,"Orthogroups.tsv")))
all_OG_long <- data.frame(rbindlist(pbmclapply(2:ncol(all_OG),function(col){
  
  # Subset...
  tmp <- data.frame(separate_rows(all_OG[,c(1,col)],colnames(all_OG)[col],sep=","))
  tmp[,2] <- gsub(" ","",tmp[,2])
  tmp$genome <- colnames(tmp)[2]
  colnames(tmp)[2] <- "gene_name_noGenome"
  tmp$gene_name <- paste0(tmp$genome,"_",tmp$gene_name_noGenome)
  return(tmp)
},mc.cores = n_cores)))

# Merge it with OG_sequences
OG_sequences2 <- merge(OG_sequences,all_OG_long[,c("Orthogroup","gene_name")],by="gene_name")
OG_sequences2$full_OF <- OG_sequences2$OF_gene

# Prepare GEA data -------------------------------------------------------
focal_datasets <- list.files("outputs/GEA_res")

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
  
  # Fetch the orthogroups as well...
  wza_merge <- merge(wza_merge,OG_sequences2[,c("Orthogroup","full_OF")],"full_OF")
  OG_counts1 <- table(wza_merge$Orthogroup)
  wza_merge <- wza_merge[wza_merge$Orthogroup %in% names(OG_counts1[OG_counts1 < orthogroup_cutoff]),]
  
  # Mark
  wza_merge$dataset <- dataset
  wza_merge$genome <- genome
  return(data.table(wza_merge[,c("dataset","Orthogroup","weiZ")]))
},mc.cores=6)
names(all_gea_res) <- focal_datasets


# Variation of method adapted to work with pairwise data... ---------------
convergence_matrix <- matrix(ncol=length(focal_datasets),nrow=length(focal_datasets))
rownames(convergence_matrix) <- colnames(convergence_matrix) <- focal_datasets
for(i in 1:nrow(convergence_matrix)){
  message(paste0(">>> Starting ",focal_datasets[i]))
  for(j in 1:ncol(convergence_matrix)){
    if(j > i){
      message(paste0(">>> Comparing with dataset ",j-i," of ",length(focal_datasets)-i))
      
      # Start by filtering our GEA result for the focal sets i and j
      all_gea_res_sub <- rbindlist(all_gea_res[c(i,j)])
      OG_subs <- unique(all_gea_res[[i]]$Orthogroup[all_gea_res[[i]]$Orthogroup %in% all_gea_res[[j]]$Orthogroup])
      all_gea_res_sub <- all_gea_res_sub[all_gea_res_sub$Orthogroup %in% OG_subs,]
      
      # Now we need to calculate empirical pvalues based on this trimmed set of wza score
      all_gea_res_sub <- all_gea_res_sub[,.(empPvals(.SD$weiZ,.SD$weiZ),.SD$Orthogroup),by=dataset]
      colnames(all_gea_res_sub) <- c("dataset","pvalue","Orthogroup")
      
      # Calculate the maximum pvalue found in each orthogroup per dataset
      OG_maxP <- all_gea_res_sub[,.(min(.SD$pvalue),length(.SD$pvalue)),by=.(Orthogroup,dataset)]
      colnames(OG_maxP) <- c("Orthogroup","dataset","maxP","Ngenes_per_species")
      
      # Define cutoffs
      cutoffs <- seq(0.001,0.05,0.001)
      
      # Dunn-Sidak correction
      OG_maxP$maxP_DS <- 1 - (1 - OG_maxP$maxP)^OG_maxP$Ngenes_per_species
      OG_maxP$maxP_DS[OG_maxP$maxP_DS > 1] <- 1
      
      # Fetch cutoff-specific distributions of pvalues
      chyper_scores <- pbmclapply(cutoffs,function(flex_cutoff){

        # Define significant orthogroups based on multiple-testing corrections
        OG_maxP$maxP_DS_signif  <- OG_maxP$maxP_DS <= flex_cutoff
        
        # Reformat for chyper
        OG_maxP1 <- OG_maxP[OG_maxP$dataset==unique(OG_maxP$dataset)[1],c("Orthogroup","maxP_DS_signif")]
        OG_maxP2 <- OG_maxP[OG_maxP$dataset==unique(OG_maxP$dataset)[2],c("Orthogroup","maxP_DS_signif")]
        OG_maxP_merge <- merge(OG_maxP1,OG_maxP2,by="Orthogroup")
        colnames(OG_maxP_merge) <- c("Orthogroup","signif1","signif2")

        
        # Calculate chyper_pairwise
        input_array <- array(0,c(nrow(OG_maxP_merge),2))
        input_array[,1] <-  as.numeric(OG_maxP_merge$signif1)
        input_array[,2] <-  as.numeric(OG_maxP_merge$signif2)
        
        # FInal matrix value
        pairwise_c_hyper(input_array)
        
        # # Fetch the number of "significant" orthogroups per species
        # OG_maxP_species_signif <- data.table(OG_maxP)[,.(sum(.SD$maxP_DS_signif),length(.SD$maxP_DS_signif)),by=dataset]
        # colnames(OG_maxP_species_signif) <- c("dataset","NSignif","N")
        # 
        # # Sum up within all OGs the number of significant pvals
        # OG_maxP_signif <- data.table(OG_maxP)[,.(sum(.SD$maxP_DS_signif),length(.SD$maxP_DS_signif)),by=Orthogroup]
        # colnames(OG_maxP_signif) <- c("Orthogroup","NSignif","N")
        # 
        # # Calculate the sum of dhyper with respect to dataset i
        # 
        # 
        # OG_maxP_signif$final_pval <- sapply(1:nrow(OG_maxP_signif),function(x) sum(dbinom(OG_maxP_signif$NSignif[x]:OG_maxP_signif$N[x],OG_maxP_signif$N[x],prob = flex_cutoff)))
        # 
        # out <- data.frame(Orthogroup=OG_maxP_signif$Orthogroup,
        #                   pvalue=OG_maxP_signif$final_pval,
        #                   fdr=p.adjust(OG_maxP_signif$final_pval,method="fdr"),
        #                   cutoff=flex_cutoff)
        # return(list(orthogroup_res=out,
        #             nsignif=length(out$fdr[out$fdr < 0.05])))
      },mc.cores=n_cores)
      
      # to_plot <- data.frame(cutoffs=cutoffs,
      #                       signif_OG=sapply(pval_dists,'[[',2))
      # ggplot(to_plot,aes(cutoffs,signif_OG))+
      #   geom_point()
      
      # FINAL DATASET SCORE
      convergence_matrix[i,j] <- mean(unlist(chyper_scores))
      
      # # FINAL OG SCORES...
      # final_OG_res <- data.frame(rbindlist(lapply(pval_dists,'[[',1)))
      # final_OG_res_scores <- sort(table(final_OG_res[final_OG_res$fdr < 0.05,"Orthogroup"]))
      # final_OG_res_scores <- rev(sort(table(final_OG_res[final_OG_res$pvalue < 0.05,"Orthogroup"])))
      
    }
  }
}

# Plot each of our matrices
plot_matrix <- function(matrix,cluster=F,phylo_order=OF_tree$tip.label){
  to_plot <- reshape2::melt(matrix)
  colnames(to_plot) <- c("Dataset1","Dataset2","Value")
  
  if(cluster){
    cluster_mat <- hclust(as.dist(matrix))
  }
  
  phylo_factor <- na.omit(dataset_meta)
  phylo_factor <- data.frame(rbindlist(lapply(phylo_order,function(x) phylo_factor[phylo_factor$genome == x,])))
  
  # Re order to match dataset...
  to_plot$phylo_order1 <- sapply(to_plot$Dataset1,function(x) which(phylo_factor$dataset == x))
  to_plot$phylo_order2 <- sapply(to_plot$Dataset2,function(x) which(phylo_factor$dataset == x))
  
  # Flip them
  for(i in 1:nrow(to_plot)){
    if(to_plot$phylo_order1[i] < to_plot$phylo_order2[i]){
      tmp <- to_plot$Dataset1[i]
      to_plot$Dataset1[i] <- to_plot$Dataset2[i]
      to_plot$Dataset2[i] <- tmp
    }
  }
  
  to_plot$Dataset1_F <- factor(to_plot$Dataset1,levels=phylo_factor$dataset)
  to_plot$Dataset2_F <- factor(to_plot$Dataset2,levels=phylo_factor$dataset)
  
  ggplot(na.omit(to_plot),aes(Dataset1_F,Dataset2_F,fill=Value))+
    geom_tile()+
    scale_fill_viridis(option="B")+
    theme(axis.text.x = element_text(angle=90,hjust=1))+
    theme(axis.title = element_blank())
}

# Get OF tree
OF_tree <- ape::read.tree("outputs/orthology/Results_210913_17_genomes_with_conifers_noAA_filter/SpeciesTree_rooted.txt")
plot_matrix(convergence_matrix,cluster = F)

# Compare against phylogenetic distance...
OF_brlens <- reshape2::melt(cophenetic.phylo(OF_tree))
OF_brlens$id <- apply(t(apply(OF_brlens[,c("Var1","Var2")],1,sort)),1,paste,collapse='-')

convergence_melt <- na.omit(reshape2::melt(convergence_matrix))
convergence_melt$genome1 <- NA
for(dataset1 in unique(convergence_melt$Var1)){
  convergence_melt[convergence_melt$Var1 == dataset1,"genome1"] <- na.omit(dataset_meta[dataset_meta$dataset == dataset1,"genome"])
}
convergence_melt$genome2 <- NA
for(dataset2 in unique(convergence_melt$Var2)){
  convergence_melt[convergence_melt$Var2 == dataset2,"genome2"] <- na.omit(dataset_meta[dataset_meta$dataset == dataset2,"genome"])
}
convergence_melt$id <- apply(t(apply(convergence_melt[,c("genome1","genome2")],1,sort)),1,paste,collapse='-')
colnames(convergence_melt)[1:3] <- c("dataset1","dataset2","Chyper")

# merge
convergence_melt_merge <- merge(convergence_melt,OF_brlens[,c("value","id")],by="id")

ggplot(convergence_melt_merge,aes(value,Chyper))+
  geom_point()+
  geom_smooth()
