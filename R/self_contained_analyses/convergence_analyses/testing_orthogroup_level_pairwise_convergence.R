# Analysis of repeatability within orthogroups among dataset pairs...
lib <- c("tidyr","ggtree","ape","VGAM","ggExtra","WeMix","pbmcapply","parallel","wCorr","data.table","ggplot2","viridis","ggridges","dplyr","patchwork","readr")
sapply(lib,library,character.only=T)

# Function Library...
calc_mapk <- function(rank_vector,rel_vector,k=length(rank_vector)){
  (1/k)*(rel_vector*log(k/rank_vector))
}

# Set variables
n_cores = 6

# Where are the blast results?
OG_dir <- "outputs/orthology/Results_210922_17_genomes_with_conifers_newPtrem_noAA_filter/Orthogroups/"
tree_dir <- "outputs/orthology/Results_210922_17_genomes_with_conifers_newPtrem_noAA_filter/Resolved_Gene_Trees/Resolved_Gene_Trees"
# OG_files <- list.files(OG_dir,pattern="Blast")

# Where is all the metadata...
vcf_genome_map <- read.table("metadata/vcf_genome_gff_210922_map.txt",fill=T)
dataset_meta <- vcf_genome_map[,c("V8","V5")]
colnames(dataset_meta) <- c("dataset","genome")


# Prepare the Orthofinder files... ----------------------------------------
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

# Prepare GEA results ------------------------------------------------------
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
},mc.cores=n_cores)
names(all_gea_res) <- focal_datasets


# Begin pairwise convergence assessment -----------------------------------
converge_matrix <- matrix(ncol=length(focal_datasets),nrow=length(focal_datasets))
rownames(converge_matrix) <- colnames(converge_matrix) <- focal_datasets

converge_corr_matrix <- converge_enrich_matrix <- converge_p_matrix <- converge_matrix

# To do this, we'll just go row by row and fill the matrix
Z_perms=200
wza_cutoffs=c(0.95,0.99)
orthogroup_cutoff <- 20 # Maximum number of genes from one species allowed in an orthogroup

if(!(file.exists(paste0("outputs/tmp_orthogroup_level_convergence_enrichment_matrices_",focal_climate,".rds")))){
  for(row in 1:nrow(converge_matrix)){
    # row_res <- pbmclapply(1:nrow(converge_matrix),function(row){
    
    # Fetch our focal dataset WZA results...
    dataset1 <- focal_datasets[row]
    message(paste0(">>> Starting ",dataset1))
    genome1 <- na.omit(dataset_meta[dataset_meta$dataset == dataset1,"genome"])
    wza1 <- readRDS(paste0("outputs/GEA_res/",dataset1,"/",focal_climate,"_WZA_TC_allgenes.rds"))
    wza1$gea_gene <- wza1$gene_id
    
    # Identify consistent cutoffs
    wza1_cutoffs <- quantile(wza1$weiZ,probs=wza_cutoffs)
    
    # Merge these with the relevant OF_code...
    OF_codes1 <- read.table(paste0("data/reference_genomes/",gff_OF_liftovers[genome1]),header=T)
    if(colnames(OF_codes1)[3]=="seqid"){
      OF_codes1$gea_gene <- paste0(OF_codes1$seqid,":",OF_codes1$start,"-",OF_codes1$end)
    }
    OF_codes1$full_OF <- paste0(genome_OF_codes[genome_OF_codes$genome == genome1,"OF_code"],"_",OF_codes1$OF_ID)
    wza1_merge <- merge(wza1[,c("gea_gene","weiZ")],OF_codes1[,c("gea_gene","full_OF")])
    wza1_merge$qseqid <- wza1_merge$full_OF
    
    # Fetch the orthogroups as well...
    wza1_merge <- merge(wza1_merge,OG_sequences2[,c("Orthogroup","full_OF")],"full_OF")
    OG_counts1 <- table(wza1_merge$Orthogroup)
    wza1_merge <- wza1_merge[wza1_merge$Orthogroup %in% names(OG_counts1[OG_counts1 < orthogroup_cutoff]),]
    
    # Now compare against the other dataset...
    for(col in 1:ncol(converge_matrix)){
      # col_res <- sapply(1:ncol(converge_matrix),function(col){
      if(row < col){ ## Ignore self-comparisons and only do both ways once...
        
        # Fetch our comparison dataset WZA results...
        dataset2 <- focal_datasets[col]
        message(paste0(">>> Comparing with dataset ",which(focal_datasets==dataset2)-row," of ",length(focal_datasets)-row))
        genome2 <- na.omit(dataset_meta[dataset_meta$dataset == dataset2,"genome"])
        wza2 <- readRDS(paste0("outputs/GEA_res/",dataset2,"/",focal_climate,"_WZA_TC_allgenes.rds"))
        wza2$gea_gene <- wza2$gene_id
        
        # Identify consistent cutoff
        wza2_cutoffs <- quantile(wza2$weiZ,probs=wza_cutoffs)
        
        # Merge these with the relevant OF_code...
        OF_codes2 <- read.table(paste0("data/reference_genomes/",gff_OF_liftovers[genome2]),header=T)
        if(colnames(OF_codes2)[3]=="seqid"){
          OF_codes2$gea_gene <- paste0(OF_codes2$seqid,":",OF_codes2$start,"-",OF_codes2$end)
        }
        OF_codes2$full_OF <- paste0(genome_OF_codes[genome_OF_codes$genome == genome2,"OF_code"],"_",OF_codes2$OF_ID)
        wza2_merge <- merge(wza2[,c("gea_gene","weiZ")],OF_codes2[,c("gea_gene","full_OF")])
        wza2_merge$sseqid <- wza2_merge$full_OF
        
        # Fetch the orthogroups as well...
        wza2_merge <- merge(wza2_merge,OG_sequences2[,c("Orthogroup","full_OF")],"full_OF")
        OG_counts2 <- table(wza2_merge$Orthogroup)
        wza2_merge <- wza2_merge[wza2_merge$Orthogroup %in% names(OG_counts2[OG_counts2 < orthogroup_cutoff]),]
        
        # We now want to highlight orthogroups containing "outlier" genes in each set...
        outlier_OG1 <- wza1_merge[wza1_merge$weiZ > wza1_cutoffs[1],"Orthogroup"]
        outlier_OG2 <- wza2_merge[wza2_merge$weiZ > wza2_cutoffs[1],"Orthogroup"]
        intersecting = Reduce(intersect,list(outlier_OG1,outlier_OG2))
        obs_enrich <- length(Reduce(intersect,list(outlier_OG1,outlier_OG2)))
        
        # First permute an orthogroup enrichment-based score...
        # There may be a bias here based on multiple OGs in each wza, but these may be of interest anyway...
        enrichment_perm <- unlist(pbmclapply(1:Z_perms,function(perm){
          set.seed(perm)
          
          rand_OG1 <- sample(wza1_merge$Orthogroup,length(outlier_OG1),replace = F)
          rand_OG2 <- sample(wza2_merge$Orthogroup,length(outlier_OG2),replace = F)
          
          length(Reduce(intersect,list(rand_OG1,rand_OG2)))
        },mc.cores = n_cores))
        
        converge_enrich_matrix[row,col] <- (obs_enrich-mean(enrichment_perm,na.rm=T))/sd(enrichment_perm,na.rm=T)
        
        # Rank-based analysis: Rank all orthogroups based on their max wza within each dataset...
        # Reorder
        wza1_merge_ordered <- wza1_merge[order(-wza1_merge$weiZ),]
        wza1_merge_ordered <- wza1_merge_ordered[!(duplicated(wza1_merge_ordered$Orthogroup)),]
        wza2_merge_ordered <- wza2_merge[order(-wza2_merge$weiZ),]
        wza2_merge_ordered <- wza2_merge_ordered[!(duplicated(wza2_merge_ordered$Orthogroup)),]
        
        # Remove any that can't be matched...
        to_keep <- table(c(wza1_merge_ordered$Orthogroup,wza2_merge_ordered$Orthogroup))
        to_keep <- names(to_keep[to_keep > 1])
        wza1_merge_ordered <- wza1_merge_ordered[wza1_merge_ordered$Orthogroup %in% to_keep,]
        wza2_merge_ordered <- wza2_merge_ordered[wza2_merge_ordered$Orthogroup %in% to_keep,]
        
        # Plot both ranks...
        both_ranks <- data.frame(Orthogroup=to_keep,
                                 rank1=match(to_keep,wza1_merge_ordered$Orthogroup),
                                 rank2=match(to_keep,wza2_merge_ordered$Orthogroup))
        both_ranks$diff <- (both_ranks$rank1-both_ranks$rank2)
        # hist(both_ranks$diff)
        
        # ggplot(both_ranks)+
        #   geom_segment(aes(x=1,xend=2,y=rank1,yend=rank2),alpha=0.2)
        # ggplot(both_ranks,aes(x=rank1,y=(diff)))+
        #   geom_point()+
        #   geom_density_2d()
        
        
        # DCG MAPK ----------------------------------------------------------------
        
        # From diffs, calculate DCG and MAPK?
        both_ranks1 <- both_ranks[order(both_ranks$rank1),]
        both_ranks1$R1_score <- (nrow(both_ranks1)-both_ranks1$rank2)/log(both_ranks1$rank1+1)
        DCG1 <- sum(abs(both_ranks1$diff/sqrt(both_ranks1$rank1+1)))
        both_ranks2 <- both_ranks[order(both_ranks$rank2),]
        DCG2 <- sum(abs(both_ranks2$diff/sqrt(both_ranks2$rank2+1)))
        DCG <- mean(abs(c(DCG1,DCG2)))
        
        # Get MAPK
        obs_mapk <- mean(sum(c(calc_mapk(both_ranks1$rank1,nrow(both_ranks1)-both_ranks1$rank2))),
                         sum(c(calc_mapk(both_ranks2$rank2,nrow(both_ranks2)-both_ranks2$rank1))))
        
        # both_ranks1$test <- test_mapk
        # ggplot(both_ranks1,aes(x=rank1,y=test))+geom_point()
        # 
        # # What the contribution from each quantile?
        # both_ranks1$rank_interval <- cut_interval(both_ranks1$rank1,n = 20)
        # prop_test <- data.frame(both_ranks1 %>% 
        #                           group_by(rank_interval) %>% 
        #                           summarise(prop=sum(test)) %>%
        #                           mutate(prop_final=prop/sum(both_ranks1$test)))
        # # Make some random ones...
        # rank_perm_sd <- unlist(pbmclapply(1:Z_perms,function(perm){
        #   set.seed(perm)
        #   
        #   perm_both_ranks <- data.frame(Orthogroup=to_keep,
        #                            rank1=sample(1:length(to_keep),length(to_keep)),
        #                            rank2=sample(1:length(to_keep),length(to_keep)))
        #   perm_both_ranks$diff <- (perm_both_ranks$rank1-perm_both_ranks$rank2)
        #   
        #   return(sd(perm_both_ranks$diff))
        # },mc.cores = n_cores))
        # 
        # # What is the rank Z score?
        # rank_Z <- (sd(both_ranks$diff)-mean(rank_perm))/sd(rank_perm)
        # converge_matrix[row,col] <- -1*rank_Z
        
        # rank_perm_DCG <- unlist(pbmclapply(1:Z_perms,function(perm){
        #   set.seed(perm)
        #   
        #   # perm_both_ranks <- data.frame(Orthogroup=to_keep,
        #   #                               rank1=sample(1:length(to_keep),length(to_keep)),
        #   #                               rank2=sample(1:length(to_keep),length(to_keep)))
        #   # perm_both_ranks$diff <- (perm_both_ranks$rank1-perm_both_ranks$rank2)
        #   
        #   # Randomise the wza scores...
        #   perm_wza1_merge_ordered <- wza1_merge[sample(1:nrow(wza1_merge)),]
        #   perm_wza1_merge_ordered <- perm_wza1_merge_ordered[!(duplicated(perm_wza1_merge_ordered$Orthogroup)),]
        #   perm_wza2_merge_ordered <- wza2_merge[sample(1:nrow(wza2_merge)),]
        #   perm_wza2_merge_ordered <- perm_wza2_merge_ordered[!(duplicated(perm_wza2_merge_ordered$Orthogroup)),]
        # 
        #   # Remove any that can't be matched...
        #   perm_to_keep <- table(c(perm_wza1_merge_ordered$Orthogroup,perm_wza2_merge_ordered$Orthogroup))
        #   perm_to_keep <- names(perm_to_keep[perm_to_keep > 1])
        #   perm_wza1_merge_ordered <- perm_wza1_merge_ordered[perm_wza1_merge_ordered$Orthogroup %in% perm_to_keep,]
        #   perm_wza2_merge_ordered <- perm_wza2_merge_ordered[perm_wza2_merge_ordered$Orthogroup %in% perm_to_keep,]
        # 
        #   # Plot both ranks...
        #   perm_both_ranks <- data.frame(Orthogroup=perm_to_keep,
        #                            rank1=match(perm_to_keep,perm_wza1_merge_ordered$Orthogroup),
        #                            rank2=match(perm_to_keep,perm_wza2_merge_ordered$Orthogroup))
        #   perm_both_ranks$diff <- (perm_both_ranks$rank1-perm_both_ranks$rank2)
        #   
        #   # From diffs, calculate DCG...
        #   perm_both_ranks1 <- perm_both_ranks[order(perm_both_ranks$rank1),]
        #   perm_DCG1 <- sum(abs(perm_both_ranks1$diff/sqrt(perm_both_ranks1$rank1+1)))
        #   perm_both_ranks2 <- perm_both_ranks[order(perm_both_ranks$rank2),]
        #   perm_DCG2 <- sum(abs(perm_both_ranks2$diff/sqrt(perm_both_ranks2$rank2+1)))
        # 
        #   return(mean(abs(c(perm_DCG1,perm_DCG2))))
        # },mc.cores = n_cores))
        
        rank_perm_MAPK <- unlist(pbmclapply(1:Z_perms,function(perm){
          set.seed(perm)
          
          # perm_both_ranks <- data.frame(Orthogroup=to_keep,
          #                               rank1=sample(1:length(to_keep),length(to_keep)),
          #                               rank2=sample(1:length(to_keep),length(to_keep)))
          # perm_both_ranks$diff <- (perm_both_ranks$rank1-perm_both_ranks$rank2)
          
          # Randomise the wza scores...
          perm_wza1_merge_ordered <- wza1_merge[sample(1:nrow(wza1_merge)),]
          perm_wza1_merge_ordered <- perm_wza1_merge_ordered[!(duplicated(perm_wza1_merge_ordered$Orthogroup)),]
          perm_wza2_merge_ordered <- wza2_merge[sample(1:nrow(wza2_merge)),]
          perm_wza2_merge_ordered <- perm_wza2_merge_ordered[!(duplicated(perm_wza2_merge_ordered$Orthogroup)),]
          
          # Remove any that can't be matched...
          perm_to_keep <- table(c(perm_wza1_merge_ordered$Orthogroup,perm_wza2_merge_ordered$Orthogroup))
          perm_to_keep <- names(perm_to_keep[perm_to_keep > 1])
          perm_wza1_merge_ordered <- perm_wza1_merge_ordered[perm_wza1_merge_ordered$Orthogroup %in% perm_to_keep,]
          perm_wza2_merge_ordered <- perm_wza2_merge_ordered[perm_wza2_merge_ordered$Orthogroup %in% perm_to_keep,]
          
          # Plot both ranks...
          perm_both_ranks <- data.frame(Orthogroup=perm_to_keep,
                                        rank1=match(perm_to_keep,perm_wza1_merge_ordered$Orthogroup),
                                        rank2=match(perm_to_keep,perm_wza2_merge_ordered$Orthogroup))
          perm_both_ranks$diff <- (perm_both_ranks$rank1-perm_both_ranks$rank2)
          
          # From diffs, calculate DCG...
          perm_both_ranks1 <- perm_both_ranks[order(perm_both_ranks$rank1),]
          perm_both_ranks2 <- perm_both_ranks[order(perm_both_ranks$rank2),]
          
          obs_mapk <- mean(sum(c(calc_mapk(perm_both_ranks1$rank1,nrow(perm_both_ranks1)-perm_both_ranks1$rank2))),
                           sum(c(calc_mapk(perm_both_ranks2$rank2,nrow(perm_both_ranks2)-perm_both_ranks2$rank1))))
          
        },mc.cores = n_cores))
        
        
        mapk_Z <- (obs_mapk - mean(rank_perm_MAPK))/sd(rank_perm_MAPK)
        converge_matrix[row,col] <- mapk_Z
        
        # Modifieid Kendalls Tau ---------------------------------------------------------------------------
        # From earlier, we can pull the ordered orthogroups
        all_OG_pairs <- pbmclapply(list(both_ranks1,both_ranks2),function(rank_tmp){
          rank_dt <- data.table(rank_tmp[,c("Orthogroup")])
          rank_dt[, `:=`(rank1 = 1L, rank2 = .I)] ## add interval columns for overlaps
          setkey(rank_dt, rank1, rank2)
          
          olaps <- foverlaps(rank_dt, rank_dt, type="within", which=TRUE)[xid != yid]
          
          all_pairs <- setDT(list(rank_dt$V1[olaps$xid],rank_dt$V1[olaps$yid]))
          # all_pairs <- paste0(rank_dt$V1[olaps$xid],"-",rank_dt$V1[olaps$yid])
          return(all_pairs)
        },mc.cores = 2)


  
  
  
  # --------------------------------------------------------------------------------------------------

      } # End if statement of row < col
    } # End col loop
  } # End row loop
} # End if file.exists

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

plot_matrix(converge_enrich_matrix)
plot_matrix(converge_matrix)

test <- na.omit(reshape2::melt(converge_matrix))
head(test[order(-test$value),])

# hist(enrichment_perm)

# # Find the orthogroup-neighbour of each outlier set
# wza1_neighbour_scores <- pbmclapply(wza1_outliers,function(gene){
#     print(gene)
#     tmp_OG <- wza1_merge[wza1_merge$full_OF==gene,"Orthogroup"]
#     tmp_gene_name <- OG_sequences2[OG_sequences2$OF_gene==gene,"gene_name"]
# 
#     # Fetch the tree
#     if(file.exists(paste0(tree_dir,"/",tmp_OG,"_tree.txt"))){
#       tmp_OG_tree <- read.tree(paste0(tree_dir,"/",tmp_OG,"_tree.txt"))
# 
#       # Filter for focal tips...
#       tips_to_keep <- tmp_OG_tree$tip.labe[tmp_OG_tree$tip.label %in%
#                                              c(tmp_gene_name,
#                                                OG_sequences2[OG_sequences2$OF_genome == genome_OF_codes[genome_OF_codes$genome==genome2,"OF_code"],"gene_name"])]
# 
#       if(length(tips_to_keep) > 1){
#         # Prune + Reformat
#         tmp_OG_dist <- reshape2::melt(as.matrix(distTips(keep.tip(tmp_OG_tree, tips_to_keep))))
#         tmp_OG_dist <- tmp_OG_dist[tmp_OG_dist$Var1 != tmp_OG_dist$Var2 &
#                                      tmp_OG_dist$Var1 == tmp_gene_name,]
#         colnames(tmp_OG_dist) <- c("focal_gene","gene_name","dist")
# 
#         # Now need to merge with the wza2 OFs and gene scores...
#         tmp_OG_dist <- merge(tmp_OG_dist,OG_sequences2[,c("gene_name","full_OF")],by="gene_name")
#         tmp_OG_dist <- merge(tmp_OG_dist,wza2_merge[,c("full_OF","weiZ")],by="full_OF")
# 
#         # And take weighted mean by distance
#         return(tmp_OG_dist[tmp_OG_dist$dist == min(tmp_OG_dist$dist),"weiZ"])
#       } else {
#         return(NA)
#       }
#     } else {
#       return(NA)
#     }
#   },mc.cores=n_cores)
# 
# # Convert to emp-p
# library(qvalue)
# test_empP <- empPvals(stat=unlist(wza1_neighbour_scores),stat0=wza2_merge$weiZ)
# test_empP <- data.frame(empP=sapply(na.omit(unlist(wza1_neighbour_scores)),function(x) length(wza2_merge$weiZ[wza2_merge$weiZ > x])/nrow(wza2_merge)))
# ggplot(test_empP,aes(empP))+
#   geom_histogram(bins=100)

#########################################################################################################
# # We also want to estimate the mean opposite wza scores for orthogroups ft. outlier genes...
# wza1_outliers <- wza1_merge[wza1_merge$weiZ > wza1_cutoffs[1],"full_OF"]
# wza2_outliers <- wza2_merge[wza2_merge$weiZ > wza2_cutoffs[1],"full_OF"]
# 
# # Do wza1 outliers first...
# wza1_opposite_scores <- pbmclapply(wza1_outliers,function(gene){
#   print(gene)
#   tmp_OG <- wza1_merge[wza1_merge$full_OF==gene,"Orthogroup"]
#   tmp_gene_name <- OG_sequences2[OG_sequences2$OF_gene==gene,"gene_name"]
#   
#   # Fetch the tree
#   if(file.exists(paste0(tree_dir,"/",tmp_OG,"_tree.txt"))){
#     tmp_OG_tree <- read.tree(paste0(tree_dir,"/",tmp_OG,"_tree.txt"))
#     
#     # Filter for focal tips...
#     tips_to_keep <- tmp_OG_tree$tip.labe[tmp_OG_tree$tip.label %in% 
#                                            c(tmp_gene_name,
#                                              OG_sequences2[OG_sequences2$OF_genome == genome_OF_codes[genome_OF_codes$genome==genome2,"OF_code"],"gene_name"])]
#     
#     if(length(tips_to_keep) > 1){
#       # Prune + Reformat
#       tmp_OG_dist <- reshape2::melt(as.matrix(distTips(keep.tip(tmp_OG_tree, tips_to_keep))))
#       tmp_OG_dist <- tmp_OG_dist[tmp_OG_dist$Var1 != tmp_OG_dist$Var2 &
#                                    tmp_OG_dist$Var1 == tmp_gene_name,]
#       colnames(tmp_OG_dist) <- c("focal_gene","gene_name","dist")
#       
#       # Now need to merge with the wza2 OFs and gene scores...
#       tmp_OG_dist <- merge(tmp_OG_dist,OG_sequences2[,c("gene_name","full_OF")],by="gene_name")
#       tmp_OG_dist <- merge(tmp_OG_dist,wza2_merge[,c("full_OF","weiZ")],by="full_OF")
#       
#       # And take weighted mean by distance
#       return(weighted.mean(tmp_OG_dist$weiZ,w = tmp_OG_dist$dist))
#     } else {
#       return(NA)
#     }
#   } else {
#     return(NA)
#   }
# },mc.cores=n_cores)
# 
# # Then wza2 outliers second...
# wza2_opposite_scores <- pbmclapply(wza2_outliers,function(gene){
#   print(gene)
#   tmp_OG <- wza2_merge[wza2_merge$full_OF==gene,"Orthogroup"]
#   tmp_gene_name <- OG_sequences2[OG_sequences2$OF_gene==gene,"gene_name"]
#   
#   # Fetch the tree
#   if(file.exists(paste0(tree_dir,"/",tmp_OG,"_tree.txt"))){
#     tmp_OG_tree <- read.tree(paste0(tree_dir,"/",tmp_OG,"_tree.txt"))
#     
#     # Filter for focal tips...
#     tips_to_keep <- tmp_OG_tree$tip.labe[tmp_OG_tree$tip.label %in% 
#                                            c(tmp_gene_name,
#                                              OG_sequences2[OG_sequences2$OF_genome == genome_OF_codes[genome_OF_codes$genome==genome1,"OF_code"],"gene_name"])]
#     
#     if(length(tips_to_keep) > 1){
#       # Prune + Reformat
#       tmp_OG_dist <- reshape2::melt(as.matrix(distTips(keep.tip(tmp_OG_tree, tips_to_keep))))
#       tmp_OG_dist <- tmp_OG_dist[tmp_OG_dist$Var1 != tmp_OG_dist$Var2 &
#                                    tmp_OG_dist$Var1 == tmp_gene_name,]
#       colnames(tmp_OG_dist) <- c("focal_gene","gene_name","dist")
#       
#       # Now need to merge with the wza2 OFs and gene scores...
#       tmp_OG_dist <- merge(tmp_OG_dist,OG_sequences2[,c("gene_name","full_OF")],by="gene_name")
#       tmp_OG_dist <- merge(tmp_OG_dist,wza1_merge[,c("full_OF","weiZ")],by="full_OF")
#       
#       # And take weighted mean by distance
#       return(weighted.mean(tmp_OG_dist$weiZ,w = tmp_OG_dist$dist))
#     } else {
#       return(NA)
#     }
#   } else {
#     return(NA)
#   }
# },mc.cores=n_cores)
#########################################################################################################       

