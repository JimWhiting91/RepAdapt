# Test of how to examine WZA outputs as manhattan plots and for clustering purposes...
lib <- c("Rfast","corrplot","EnvStats","zoo","data.table","tidyverse","ggplot2","gtools","pbmcapply","regioneR","cowplot","pracma")
lapply(lib,library,character.only=T)

# Fetch all of our output dir
result_dir <- list.files("outputs/GEA_res")

# For now we'll look at Amaranthus
species_res <- "Amaranthus_tuberculatus_Wright_Individual"
gea_res <- list.files(paste0("outputs/GEA_res/",species_res,"/"),pattern = "WZA_empiricalp_full")

# Set wingspan for overlap
wingspan_size = 0.5e6

# Fetch the fai
fai <- read.table("data/reference_genomes/A.tuberculatus_reference/Amaranthus_tuberculatus.fasta.fai")

# Run over AMT for now
all_res <- data.frame(rbindlist(pbmclapply(gea_res,function(res){
  tmp_res <- data.frame(fread(paste0("outputs/GEA_res/",species_res,"/",res)))
  
  # Reconstruct the chromosome and start/end positions
  chr <- sapply(strsplit(tmp_res$gene_id,":"),'[[',1)
  pos <- sapply(strsplit(tmp_res$gene_id,":"),'[[',2)
  start <- as.integer(sapply(strsplit(pos,"-"),'[[',1))
  end <- as.integer(sapply(strsplit(pos,"-"),'[[',2))
  tmp_res$chr <- chr
  tmp_res$start <- start
  tmp_res$end <- end
  tmp_res$mid <- rowMeans(tmp_res[,c("start","end")])
  tmp_res$chr_F <- factor(tmp_res$chr,levels=mixedsort(unique(tmp_res$chr)))
  
  return(tmp_res)
},mc.cores = 4)))

# Get colours
scaf_cols <- data.frame(chr_F=mixedsort(unique(all_res$chr)),
                        cols=rep(c("black","grey50"),length(unique(all_res$chr))/2))

# Get our cutoff thresholds
cutoff_thresholds <- data.frame(climate_var=unique(all_res$climate_var))
cutoff_thresholds$alpha05 <- NA
cutoff_thresholds$alpha01 <- NA
for(var in cutoff_thresholds$climate_var){
  cutoff_thresholds[cutoff_thresholds$climate_var == var,"alpha05"] <- quantile(all_res[all_res$climate_var == var,"weiZ"],0.95)
  cutoff_thresholds[cutoff_thresholds$climate_var == var,"alpha01"] <- quantile(all_res[all_res$climate_var == var,"weiZ"],0.99)
}

# Try a test plot
#genome_weiZ <- 
ggplot(all_res[all_res$climate_var %in% c("isothermality","temp_seasonality"),],aes(mid,weiZ,colour=chr_F))+
  geom_point(show.legend = F)+
  #geom_line()+
  facet_grid(climate_var ~ chr_F, scales = "free", space='free')+
  theme(panel.spacing = unit(0, "lines"))+
  scale_colour_manual(breaks=scaf_cols$chr_F,
                      values=scaf_cols$cols)


# Regioner approach --------------------------------------------------------
# Set up our GRanges for all genes
just_genes <- unique(all_res[,c("chr","start","end","gene_id")])
all_gene_grange <- GRanges(seqnames = just_genes$chr,
                           ranges = IRanges(start = just_genes$start,
                                            end = just_genes$end,
                                            names = just_genes$gene_id))

# Set up genome
genome_fai <- data.frame(chr=fai[,1],
                         start=1,
                         end=fai[,2])

# Now run over individual climate_var
climate_var_clustering <- lapply(unique(all_res$climate_var),function(climate_var){
  
  print(paste0("STARTING ",climate_var))
  
  # Make our focal outlier set
  tmp_res <- all_res[all_res$climate_var == climate_var,]
  focal_genes <- tmp_res[tmp_res$weiZ >= quantile(tmp_res$weiZ,0.99),c("chr","start","end","gene_id")]
  focal_genes_range <- GRanges(seqnames = focal_genes$chr,
                               ranges = IRanges(start = focal_genes$start,
                                                end = focal_genes$end,
                                                names = focal_genes$gene_id))
  
  # Remove our focal genes from the all genes universe
  neutral_genes <- tmp_res[tmp_res$weiZ < quantile(tmp_res$weiZ,0.99),c("chr","start","end","gene_id")]
  neutral_genes_range <- GRanges(seqnames = neutral_genes$chr,
                                 ranges = IRanges(start = neutral_genes$start,
                                                  end = neutral_genes$end,
                                                  names = neutral_genes$gene_id))
  
  # Now we want to calculate the distance between genes on the same chromosome...
  focal_nearest <- distanceToNearest(focal_genes_range,ignore.strand=T)
  focal_nearest_res <- as.data.frame(focal_nearest)
  
  # Also count the proportion of non-overlapping genes that overlap following wingspan adjustments
  focal_genes_range2 <- joinRegions(focal_genes_range)
  
  # Now extend wingspans...
  focal_genes_range2 <- extendRegions(focal_genes_range2,
                                      extend.start = wingspan_size,
                                      extend.end = wingspan_size)
  joined_focal_regions <- joinRegions(focal_genes_range2, min.dist=1)
  
  # Whats the count of overlapping regions...
  overlap_prop <- length(joined_focal_regions)/length(focal_genes_range2)
  
  # And now we want to collect a random permutation of distances/overlap among genes...
  perms=1000
  perm_distances_overlap <- pbmclapply(1:perms,function(iter){
    
    # Resample genes from the all genes universe
    perm_regions <- resampleRegions(focal_genes,genome=genome_fai,universe=neutral_genes_range, per.chromosome=FALSE)
    
    # Fetch the nearest distance
    perm_dist <- as.data.frame(distanceToNearest(perm_regions,ignore.strand=T))
    perm_dist$iter <- iter
    
    # And get overlap
    # Now extend wingspans...
    perm_regions2 <- extendRegions(perm_regions,
                                   extend.start = wingspan_size,
                                   extend.end = wingspan_size)
    joined_perm_regions <- joinRegions(perm_regions2, min.dist=1)
    
    # Whats the count of overlapping regions...
    perm_overlap_prop <- length(joined_perm_regions)/length(perm_regions2)
    
    return(list(perm_dist,perm_overlap_prop))
  },mc.cores=6)
  
  # Separate back out
  perm_distances <- lapply(perm_distances_overlap,'[[',1)
  perm_overlap <- sapply(perm_distances_overlap,'[[',2)
  
  # From our perms, we are then interested in summarising each permuted distributions
  all_medians <- sapply(perm_distances,function(x) median(x$distance))
  obs_Z <- (median(focal_nearest_res$distance) - median(all_medians))/sd(all_medians)
  obs_p <- length(all_medians[all_medians <= median(focal_nearest_res$distance)])/length(all_medians)
  
  # Return these results for the climate variable
  out_res <- data.frame(climate_var = climate_var,
                        median_dist = median(focal_nearest_res$distance),
                        distance_obs_Z = obs_Z,
                        distance_obs_p = obs_p,
                        wingspan_overlap = overlap_prop,
                        wingspan_obs_Z = (overlap_prop - mean(perm_overlap))/sd(perm_overlap),
                        wingspan_obs_p = length(perm_overlap[perm_overlap <= overlap_prop])/length(perm_overlap))
  
  return(list(out_res,all_medians))
})

# Make a summary figure...
all_observed <- data.frame(rbindlist(lapply(climate_var_clustering,'[[',1)))
all_perms <- data.frame(perm_dist=unlist(lapply(climate_var_clustering,'[[',2)),
                        climate_var=rep(unique(all_res$climate_var),each=length(climate_var_clustering[[1]][[2]])))
all_perms$obs_dist <- rep(all_observed$median_dist,each=nrow(all_perms)/length(unique(all_perms$climate_var)))

# Figure summarising the Z-scores of each clustering approach...
ggplot(all_observed,aes(abs(distance_obs_Z),abs(wingspan_obs_Z)))+
  geom_point()+
  geom_abline()+
  labs(y="Wingspan Overlap",x="Median Distance to Neighbour")

# Figure summary
perm_summary_fig <- ggplot(all_perms[all_perms$climate_var %in% unique(all_perms$climate_var)[1:5],],aes(perm_dist))+
  geom_histogram(bins=50,fill="gray50")+
  geom_vline(aes(xintercept = obs_dist),colour="red2")+
  # geom_vline(xintercept = median(all_medians),colour="black",linetype="dashed")+
  theme_minimal()+
  facet_wrap(~climate_var,ncol=1,strip.position = "right")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=14))+
  labs(x="Distance to Nearest Neighbour (Mb)")+
  scale_x_continuous(labels = function(x) x/1000000)

# Combine the summary with the manhattans
plot_grid(genome_weiZ,perm_summary_fig,
          rel_widths = c(3,1),align="h",ncol=2,nrow=1,axis="tblr")

# Just plot all of the climate_var distances...
all_perm_res <- ggplot(all_perms,aes(perm_dist))+
  geom_histogram(bins=50,fill="gray50")+
  geom_vline(aes(xintercept = obs_dist),colour="red2")+
  # geom_vline(xintercept = median(all_medians),colour="black",linetype="dashed")+
  theme_minimal()+
  facet_wrap(~climate_var,ncol=5,strip.position = "right")+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=14))+
  labs(x="Distance to Nearest Neighbour (Mb)")+
  scale_x_continuous(labels = function(x) x/1000000)

# And all the Z scores
Zscore_plot <- ggplot(all_observed,aes(y=climate_var,x=obs_Z))+
  geom_bar(stat="identity")
plot_grid(all_perm_res,Zscore_plot,
          ncol=2,rel_widths=c(2,1))

# What associates with NN Z-scores
climate_data <- read.table(paste0("outputs/GEA_res/",species_res,"/climate_cline.tsv"),header=T)

# What are CV for each
clim_CV <- apply(climate_data[,c(unique(all_res$climate_var))],2,cv)
clim_PCA <- prcomp(climate_data[,c(unique(all_res$climate_var))],scale=T,center = T)
corrplot(cor(climate_data[,c(unique(all_res$climate_var))]))

# Associated with 'peak height'
cor.test(cutoff_thresholds$alpha01,all_observed$Obs_Z)


# Per-Chromosome Clustering... --------------------------------------------
# Repeat of overlap-based clustering, but derived within chromosomes and a Z score given for each
# Now run over individual climate_var

chr_climate_var_clustering <- lapply(unique(all_res$climate_var),function(climate_var){
  
  print(paste0("STARTING ",climate_var))
  
  # Make our focal outlier set
  tmp_res <- all_res[all_res$climate_var == climate_var,]
  focal_genes <- tmp_res[tmp_res$weiZ >= quantile(tmp_res$weiZ,0.99),c("chr","start","end","gene_id")]
  focal_genes_range <- GRanges(seqnames = focal_genes$chr,
                               ranges = IRanges(start = focal_genes$start,
                                                end = focal_genes$end,
                                                names = focal_genes$gene_id))
  
  # Remove our focal genes from the all genes universe
  neutral_genes <- tmp_res[tmp_res$weiZ < quantile(tmp_res$weiZ,0.99),c("chr","start","end","gene_id")]
  neutral_genes_range <- GRanges(seqnames = neutral_genes$chr,
                                 ranges = IRanges(start = neutral_genes$start,
                                                  end = neutral_genes$end,
                                                  names = neutral_genes$gene_id))
  
  # Now we want to calculate the distance between genes on the same chromosome...
  focal_nearest <- distanceToNearest(focal_genes_range,ignore.strand=T)
  focal_nearest_res <- as.data.frame(focal_nearest)
  
  # Also count the proportion of non-overlapping genes that overlap following wingspan adjustments
  focal_genes_range2 <- joinRegions(focal_genes_range)
  
  # Now extend wingspans...
  focal_genes_range2 <- extendRegions(focal_genes_range2,
                                      extend.start = wingspan_size,
                                      extend.end = wingspan_size)
  joined_focal_regions <- joinRegions(focal_genes_range2, min.dist=1)
  
  # Now divide these up into chromosomal segments...
  chr_focal_genes_range2 <- split(focal_genes_range2, mixedsort(seqnames(focal_genes_range2)))
  chr_joined_focal_regions <- split(joined_focal_regions, mixedsort(seqnames(joined_focal_regions)))
  
  # Whats the count of overlapping regions...
  chr_overlap_prop <- matrix(ncol=length(chr_joined_focal_regions),nrow=1)
  colnames(chr_overlap_prop) <- names(chr_joined_focal_regions)
  for(i in 1:ncol(chr_overlap_prop)){
    chr_overlap_prop[1,i] <- length(chr_joined_focal_regions[[i]])/length(chr_focal_genes_range2[[i]])
  }
  
  # And now we want to collect a random permutation of distances/overlap among genes...
  perms=1000
  perm_chr_overlap <- data.frame(rbindlist(pbmclapply(1:perms,function(iter){
    
    # Resample genes from the all genes universe
    perm_regions <- resampleRegions(focal_genes,genome=genome_fai,universe=neutral_genes_range, per.chromosome=FALSE)
    
    # And get overlap
    # Now extend wingspans...
    perm_genes_regions2 <- extendRegions(perm_regions,
                                         extend.start = wingspan_size,
                                         extend.end = wingspan_size)
    joined_perm_regions <- joinRegions(perm_genes_regions2, min.dist=1)
    
    # Sort
    perm_genes_regions2 <- sort(perm_genes_regions2)
    joined_perm_regions <- sort(joined_perm_regions)
    
    # Now divide these up into chromosomal segments...
    chr_perm_genes_regions2 <- split(perm_genes_regions2, mixedsort(seqnames(perm_genes_regions2)))
    chr_joined_perm_regions <- split(joined_perm_regions, mixedsort(seqnames(joined_perm_regions)))
    
    # Whats the count of overlapping regions...
    perm_chr_overlap_prop <- matrix(ncol=length(chr_joined_perm_regions),nrow=1)
    colnames(perm_chr_overlap_prop) <- names(chr_joined_perm_regions)
    for(i in 1:ncol(perm_chr_overlap_prop)){
      perm_chr_overlap_prop[1,i] <- length(chr_joined_perm_regions[[i]])/length(chr_perm_genes_regions2[[i]])
    }
    
    return(data.frame(perm_chr_overlap_prop))
  },mc.cores=6)))
  
  # Can now extract Z-scores...
  out <- data.frame(t(chr_overlap_prop))
  colnames(out) <- "Obs_overlap"
  
  # Fetch Z-scores
  out$Z_overlap <- NA
  for(i in 1:nrow(out)){
    out$Z_overlap[i] <- (out$Obs_overlap[i] - mean(perm_chr_overlap[,i]))/sd(perm_chr_overlap[,i])
  }
  
  # Add identifiers...
  out$chr <- rownames(out)
  out$climate_var <- climate_var
  
  return(out)
})

# Combine all these together
all_chr_clusters <- data.frame(rbindlist(chr_climate_var_clustering))

# What's the variance for each climate var...
all_chr_clusters %>% group_by(climate_var) %>% summarise(chr_cluster_var=var(Z_overlap))

# Also plot...
all_chr_clusters$chr_F <- factor(all_chr_clusters$chr,levels=rev(mixedsort(unique(all_chr_clusters$chr))))
ggplot(all_chr_clusters,aes(x=Z_overlap,y=chr_F))+
  geom_bar(stat = "identity")+
  facet_wrap(~climate_var)+
  geom_vline(xintercept = -3,colour="red2")

# And calculate the variance...
chr_cluster_variance <- all_chr_clusters %>% group_by(climate_var) %>% summarise(cluster_var = var(Z_overlap))
chr_cluster_variance <- chr_cluster_variance[order(-chr_cluster_variance$cluster_var),]
all_chr_clusters$climate_var_F <- factor(all_chr_clusters$climate_var,levels=chr_cluster_variance$climate_var)
ggplot(all_chr_clusters,aes(x=Z_overlap,y=chr_F))+
  geom_bar(stat = "identity")+
  facet_wrap(~climate_var_F)+
  geom_vline(xintercept = -3,colour="red2")

# Comparison of clustering among outlier sets across environmental --------

# First fetch our environmental variables and build correlation matrix
env_variables <- read.table("outputs/GEA_res/Amaranthus_tuberculatus_Wright_Individual/climate_cline.tsv",header=T)
env_variables <- env_variables[,unique(all_res$climate_var)]
env_corr_mat <- cor(env_variables,method = "spearman")

# Build a list of all focal outlier sets...
focal_genes_list <- lapply(unique(all_res$climate_var),function(climate_var){
  
  # Make our focal outlier set
  tmp_res <- all_res[all_res$climate_var == climate_var,]
  focal_genes <- tmp_res[tmp_res$weiZ >= quantile(tmp_res$weiZ,0.99),c("chr","start","end","gene_id")]
  focal_genes_range <- GRanges(seqnames = focal_genes$chr,
                               ranges = IRanges(start = focal_genes$start,
                                                end = focal_genes$end,
                                                names = focal_genes$gene_id))
  
  # Add extensions
  focal_genes_range2 <- joinRegions(focal_genes_range)
  focal_genes_range2 <- extendRegions(focal_genes_range2,
                                      extend.start = wingspan_size,
                                      extend.end = wingspan_size)
  return(focal_genes_range2)
})

# Convert to GRangesList
focal_genes_list <- GRangesList(focal_genes_list[1:length(focal_genes_list)])

# Estimate overlap among all pairwise...
overlap_matrix <- matrix(ncol=ncol(env_corr_mat),nrow=nrow(env_corr_mat),NA)
rownames(overlap_matrix) <- rownames(env_corr_mat)
colnames(overlap_matrix) <- colnames(env_corr_mat)

# Add wingspans onto random genes
all_gene_grange_wings <- extendRegions(all_gene_grange,
                                       extend.start = wingspan_size,
                                       extend.end = wingspan_size)

# Populate with the overlap...
for(i in 1:nrow(overlap_matrix)){
  print(paste0("Filling row ",rownames(overlap_matrix)[i]))
  for(j in 1:ncol(overlap_matrix)){
    
    # Only fill one half...
    if(j > i){
      print(paste0("Running ",j," of ",ncol(overlap_matrix)))
      
      # Count the overlap among the ith and jth elements
      overlap_Z <- permTest(A=focal_genes_list[[i]],
                            B=focal_genes_list[[j]],
                            ntimes = 1000,
                            alternative = "greater",
                            randomize.function=resampleRegions,
                            universe=all_gene_grange_wings,
                            evaluate.function = numOverlaps,
                            force.parallel = T)
      
      # populate the matrix
      overlap_matrix[i,j] <- overlap_Z$numOverlaps$zscore
    }
  }
}

# Plot them side by side...
overlap_melt <- reshape2::melt(overlap_matrix)
corr_melt <- reshape2::melt(env_corr_mat)
overlap_melt$env_corr <- corr_melt$value

# Remove lower tri
overlap_melt <- na.omit(overlap_melt)

# Remove self-overlap...
overlap_melt <- overlap_melt[overlap_melt$Var1 != overlap_melt$Var2,]
colnames(overlap_melt) <- c("Var1","Var2","Overlap_Z","Env_Corr")
overlap_melt$Var1_F <- factor(overlap_melt$Var1,levels=unique(all_res$climate_var))
overlap_melt$Var2_F <- factor(overlap_melt$Var2,levels=unique(all_res$climate_var))

# Plot both
plot_grid(
ggplot(overlap_melt,aes(Var1_F,Var2_F,fill=Overlap_Z))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme(axis.text.x = element_text(angle=45,hjust=1)),
ggplot(overlap_melt,aes(Var1_F,Var2_F,fill=abs(Env_Corr)))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme(axis.text.x = element_text(angle=45,hjust=1)),
ncol=2)

# Also plot them against one another
ggplot(overlap_melt,aes(abs(Env_Corr),Overlap_Z))+
  geom_point()

# # Melt again...
# overlap_melt2 <- reshape2::melt(overlap_melt)
# 
# # Plot
# overlap_melt2$Var1_F <- factor(overlap_melt2$Var1,levels=unique(all_res$climate_var))
# overlap_melt2$Var2_F <- factor(overlap_melt2$Var2,levels=unique(all_res$climate_var))
# ggplot(overlap_melt2,aes(Var1_F,Var2_F,fill=value))+
#   geom_tile()+
#   facet_wrap(~variable,ncol=2)
