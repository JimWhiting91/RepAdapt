# Pairwise orthology based on OF2 Blast results...
lib <- c("qvalue","wrswoR","ggExtra","pbmcapply","parallel","data.table","ggplot2","viridis","ggridges","dplyr","patchwork","readr")
sapply(lib,library,character.only=T)

# Where are the blast results?
blast_dir <- "outputs/orthology/Results_210913_17_genomes_with_conifers_noAA_filter/WorkingDirectory/"
blast_files <- list.files(blast_dir,pattern="Blast")

# Where is all the metadata...
vcf_genome_map <- read.table("metadata/vcf_genome_gff_210922_map.txt")
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

# Let's converge... -------------------------------------------------------
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

converge_matrix <- matrix(ncol=length(focal_datasets),nrow=length(focal_datasets))
rownames(converge_matrix) <- colnames(converge_matrix) <- focal_datasets

converge_corr_matrix <- converge_enrich_matrix <- converge_p_matrix <- converge_matrix

# To do this, we'll just go row by row and fill the matrix
Z_perms=100
wza_cutoffs=c(0.95,0.99)

# # Make our output matrices...
# converge_enrich_matrix_list <- lapply(1:length(wza_cutoffs),function(i) matrix(ncol=length(focal_datasets),nrow=length(focal_datasets)))
# names(converge_enrich_matrix_list) <- paste0("Cutoff = ",wza_cutoffs)


if(!(file.exists(paste0("outputs/tmp_blast-aware_convergence_enrichment_matrices_",focal_climate,".rds")))){
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
        
        # Read in both blast_results...
        blast_res1 <- suppressMessages(
          data.frame(read_table(paste0(blast_dir,"/Blast",genome_OF_codes[genome_OF_codes$genome == genome1,"OF_code"],"_",genome_OF_codes[genome_OF_codes$genome == genome2,"OF_code"],".txt.gz"),progress = F,col_names = F))
        )
        colnames(blast_res1) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","start","ssend","evalue","bitscore")
        blast_res2 <- suppressMessages(
          data.frame(read_table(paste0(blast_dir,"/Blast",genome_OF_codes[genome_OF_codes$genome == genome2,"OF_code"],"_",genome_OF_codes[genome_OF_codes$genome == genome1,"OF_code"],".txt.gz"),progress = F,col_names = F))
        )
        colnames(blast_res2) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","start","ssend","evalue","bitscore")
        
        # Reciprocally filter both of these
        blast_res1$qs_ortho <- paste0(blast_res1$qseqid,"-",blast_res1$sseqid)
        blast_res1$sq_ortho <- paste0(blast_res1$sseqid,"-",blast_res1$qseqid)
        blast_res2$qs_ortho <- paste0(blast_res2$qseqid,"-",blast_res2$sseqid)
        blast_res2$sq_ortho <- paste0(blast_res2$sseqid,"-",blast_res2$qseqid)
        blast_res1_filt <- blast_res1[blast_res1$qs_ortho %in% blast_res2$sq_ortho,]
        blast_res2_filt <- blast_res2[blast_res2$qs_ortho %in% blast_res1$sq_ortho,]
        
        # Add dataset IDs
        blast_res1_filt$dataset <- dataset1
        blast_res2_filt$dataset <- dataset2
        
        # Flip the qseq-sseq for res2 and merge all the results into a single file...
        tmp <- blast_res2_filt$qseqid
        blast_res2_filt$qseqid <- blast_res2_filt$sseqid
        blast_res2_filt$sseqid <- tmp
        
        blast_res_filt <- rbind(blast_res1_filt,blast_res2_filt)
        
        # # Transform evalues
        # blast_res$evalue_trans <- -log10(blast_res$evalue + min(blast_res[blast_res$evalue != 0,"evalue"]))
        
        # Attempt at "correlation of repeatability" approach ----------------------
        
        # And we want the ranking of every gene in both seq columns...
        blast_res_clean <- merge(blast_res_filt,wza1_merge[,c("qseqid","weiZ")],"qseqid")
        colnames(blast_res_clean)[ncol(blast_res_clean)] <- "q_wza"
        blast_res_clean <- merge(blast_res_clean,wza2_merge[,c("sseqid","weiZ")],"sseqid")
        colnames(blast_res_clean)[ncol(blast_res_clean)] <- "s_wza"
        blast_res_clean <- blast_res_clean[,c("qseqid","sseqid","bitscore","q_wza","s_wza","dataset")]
        
        # # Now perform a weighted correlation...
        # # test_corr <- weightedCorr(blast_res_clean$q_wza,blast_res_clean$s_wza,method = "Spearman",weights = blast_res_clean$bitscore)
        # upper_half <- blast_res_clean[blast_res_clean$q_wza > 0 & blast_res_clean$s_wza > 0,]
        # upper_half_corrs <- sapply(c(0.5,0.9,0.95,0.99),function(prob){
        #   
        #   # Fetch new wza cutoffs
        #   wza1_cutoff_tmp <- quantile(unique(blast_res_clean[,c("qseqid","q_wza")])[,"q_wza"],probs=prob)
        #   wza2_cutoff_tmp <- quantile(unique(blast_res_clean[,c("sseqid","s_wza")])[,"s_wza"],probs=prob)
        #   
        #   # Weighted correlation on this subset of the data
        #   sub_blast_res <- blast_res_clean[blast_res_clean$q_wza > wza1_cutoff_tmp &
        #                                      blast_res_clean$s_wza > wza2_cutoff_tmp,]
        #   
        #   ggplot(sub_blast_res,aes(q_wza,s_wza,colour=bitscore))+geom_point()
        #   
        #   # Get corr
        #   weightedCorr(sub_blast_res$q_wza,sub_blast_res$s_wza,method = "Spearman",weights = sub_blast_res$bitscore)
        # })
        
        # Build a weighted mixed-effect model...
        # m1 <- mix(s_wza ~ q_wza + (sseqid|qseqid), data=upper_half, weights="evalue_trans")
        # m1 <- lmer(q_wza ~ s_wza*evalue_trans + (sseqid|qseqid), data=upper_half)
        # plot(allEffects(m1))
        #         
        #         test <- ggplot(blast_res_clean,aes(q_wza,s_wza))+geom_point()+
        #           geom_vline(xintercept = wza1_cutoff,colour="blue2")+
        #           geom_hline(yintercept = wza2_cutoff,colour="blue2")+
        #           geom_vline(xintercept = 0,colour="red2")+
        #           geom_hline(yintercept = 0,colour="red2")
        #         
        #         ggMarginal(test, type="histogram")
        #         
        #         return(upper_corr)
        #       }
        #     })
        # return(col_res)
        
        # Attempt at "permuted enrichment of repeatability" approach ----------------------
        # Do this once to prepare for permutations
        wza_outliers <- lapply(wza_cutoffs,function(cutoff){
          wza1_outliers <- wza1_merge[wza1_merge$weiZ > quantile(wza1_merge$weiZ,cutoff),]
          wza2_outliers <- wza2_merge[wza2_merge$weiZ > quantile(wza2_merge$weiZ,cutoff),]
          return(list(wza1=wza1_outliers,wza2=wza2_outliers))
        })
        
        # # Just look for 'any' match...
        # wza1_orthologs <- blast_res_clean[blast_res_clean$qseqid %in% wza1_outliers$qseqid,]
        # 
        # # Add some bit groupings
        # wza1_orthologs$discrete_eval <- cut_interval(wza1_orthologs$evalue_trans,n=10)
        # 
        # ggplot(wza1_orthologs,aes(q_wza,s_wza))+
        #   geom_point()
        # ggplot(wza1_orthologs,aes(x=s_wza,y=discrete_eval))+
        #   geom_violin(draw_quantiles = 0.5)
        
        # # It might be quicker to just expand out wza1_outliers once rather than repeatedly sample within the loops...
        # blast_res_expanded <- data.frame(rbindlist(lapply(unique(wza1_outliers$qseqid),function(qseqid){
        #   tmp <- blast_res_clean[blast_res_clean$qseqid==qseqid,]
        #   tmp[rep(seq_len(nrow(tmp)), times = tmp$bitscore), ]
        # })))
        # 
        # perm_res <- mclapply(1:Z_perms,function(iter){
        #   set.seed(iter)
        #   
        #   # Take a sample from the expanded blast_res_clean
        #   sample_tmp <- blast_res_expanded[sample(1:nrow(blast_res_expanded)),]
        #   sample_tmp <- sample_tmp[!(duplicated(sample_tmp$qseqid)),]
        #   tmp_ortho_fetch <- sample_tmp$sseqid
        #   
        #   # # How many of these are similarly in the top 1% - We expect 1%...
        #   # length(tmp_ortho_fetch[tmp_ortho_fetch %in% wza2_outliers$sseqid])/nrow(wza2_outliers)
        #   
        #   # How shifted is the distribution between these and without
        #   chosen_orthos <- wza2_merge[wza2_merge$sseqid %in% tmp_ortho_fetch,]
        #   other_orthos <- wza2_merge[!(wza2_merge$sseqid %in% tmp_ortho_fetch) &
        #                                wza2_merge$sseqid %in% unique(blast_res_clean$sseqid),]
        #   ortho_Z <- (mean(chosen_orthos$weiZ)-mean(other_orthos$weiZ))/sqrt(var(chosen_orthos$weiZ) / nrow(chosen_orthos) + var(other_orthos$weiZ) / nrow(other_orthos))
        #   return(ortho_Z)
        # },mc.cores=6)
        
        # # Permutations with C-based resampling...
        # perm_res <- data.frame(rbindlist(pbmclapply(1:Z_perms,function(iter){
        #   set.seed(iter)
        #   
        #   # Faster weighted sampling...
        #   sample_tmp <- blast_res_clean[sample_int_crank(nrow(blast_res_clean),size=nrow(blast_res_clean),prob = blast_res_clean$bitscore),]
        #   
        #   # Run our resampling over all cutoffs
        #   cutoff_enrichments <- sapply(1:length(wza_cutoffs),function(cutoff){
        #     
        #     # Draw and condense down our wza-outlier samples
        #     wza1_sampled <- sample_tmp[sample_tmp$qseqid %in% wza_outliers[[cutoff]]$wza1$full_OF &
        #                                  sample_tmp$dataset == dataset1,]
        #     wza1_sampled <- wza1_sampled[!(duplicated(wza1_sampled$qseqid)),]
        #     wza2_sampled <- sample_tmp[sample_tmp$sseqid %in% wza_outliers[[cutoff]]$wza2$full_OF &
        #                                  sample_tmp$dataset == dataset2,]
        #     wza2_sampled <- wza2_sampled[!(duplicated(wza2_sampled$sseqid)),]
        #     
        #     # Draw and condense down our random samples
        #     random1_sampled <- sample_tmp[!(sample_tmp$qseqid %in% wza_outliers[[cutoff]]$wza1$full_OF) &
        #                                     sample_tmp$dataset == dataset1,]
        #     random1_sampled <- random1_sampled[!(duplicated(random1_sampled$qseqid)),]
        #     random1_sampled <- random1_sampled[1:nrow(wza1_sampled),]
        #     
        #     random2_sampled <- sample_tmp[!(sample_tmp$sseqid %in% wza_outliers[[cutoff]]$wza2$full_OF) &
        #                                     sample_tmp$dataset == dataset2,]
        #     random2_sampled <- random2_sampled[!(duplicated(random2_sampled$sseqid)),]
        #     random2_sampled <- random2_sampled[1:nrow(wza2_sampled),]
        #     
        #     # # Compare the overlap of each set - unique? or not?
        #     # wza_overlap <- mean(c(length(which(unique(wza1_sampled$sseqid) %in% unique(wza2_sampled$sseqid))), length(which(unique(wza1_sampled$qseqid) %in% unique(wza2_sampled$qseqid)))))
        #     # random_overlap <- mean(c(length(which(unique(random1_sampled$sseqid) %in% unique(random2_sampled$sseqid))), length(which(unique(random1_sampled$qseqid) %in% unique(random2_sampled$qseqid)))))
        #     
        #     # Compare the overlap of each set - unique? or not?
        #     enrich1 <- (length(which(unique(wza1_sampled$sseqid) %in% unique(wza2_sampled$sseqid)))+1)/(length(which(unique(random1_sampled$sseqid) %in% unique(wza2_sampled$sseqid)))+1)
        #     enrich2 <- (length(which(unique(wza2_sampled$qseqid) %in% unique(wza1_sampled$qseqid)))+1)/(length(which(unique(random2_sampled$qseqid) %in% unique(wza1_sampled$qseqid)))+1)
        #     
        #     # return(mean(c(enrich1,enrich2)))
        #     # Return the geometric mean of 1 and 2
        #     return(exp(mean(log(c(enrich1,enrich2)))))
        #   })
        #   
        #   data.frame(cutoff=wza_cutoffs,
        #              enrichment=cutoff_enrichments)
        # },mc.cores=6)))
        # 
        # # Visualise
        # ggplot(perm_res,aes(x=enrichment))+
        #   geom_histogram()+
        #   facet_wrap(~cutoff)
        # 
        # # Summarise
        # for(i in 1:length(wza_cutoffs)){
        #   converge_enrich_matrix_list[[i]][row,col] <- mean(perm_res[perm_res$cutoff == wza_cutoffs[i],"enrichment"],na.rm=T)
        # }
        # 
        # converge_enrich_matrix[row,col] <- mean(unlist(perm_res))
        
        
        # Approach of weighted mean per orthogroup... -----------------------------
        # Calculate the weighted wza for all orthologs in outliers/non-outliers
        wza1_means <- data.frame(blast_res_clean[blast_res_clean$dataset==dataset1,] %>% group_by(qseqid) %>% summarise(blast_wza=weighted.mean(s_wza,w=bitscore)))
        wza2_means <- data.frame(blast_res_clean[blast_res_clean$dataset==dataset2,] %>% group_by(sseqid) %>% summarise(blast_wza=weighted.mean(q_wza,w=bitscore)))
        
        # Merge with outlier status...
        wza1_means$outlier <- "Neutral"
        wza2_means$outlier <- "Neutral"
        # for(i in 1:length(wza_cutoffs)){
        #   wza1_means[wza1_means$qseqid %in% wza_outliers[[i]]$wza1$full_OF,"outlier"] <- wza_cutoffs[i]
        #   wza2_means[wza2_means$sseqid %in% wza_outliers[[i]]$wza2$full_OF,"outlier"] <- wza_cutoffs[i]
        # }
        
        # Define outliers as top 10%...
        wza1_means[wza1_means$qseqid %in% wza_outliers[[2]]$wza1$full_OF,"outlier"] <- "Top5"
        wza2_means[wza2_means$sseqid %in% wza_outliers[[2]]$wza2$full_OF,"outlier"] <- "Top5"
        
        # Calculate empirical p-vals for outliers...
        wza1_means_outliers <- wza1_means[!(wza1_means$outlier %in% c("Neutral")),]
        wza2_means_outliers <- wza2_means[!(wza2_means$outlier %in% c("Neutral")),]
        
        wza1_means_outliers$emp_p <- empPvals(wza1_means_outliers$blast_wza,stat0=wza1_means[wza1_means$outlier=="Neutral","blast_wza"])
        wza2_means_outliers$emp_p <- empPvals(wza2_means_outliers$blast_wza,stat0=wza2_means[wza2_means$outlier=="Neutral","blast_wza"])
        
        # Calc FDR
        wza1_means_outliers$fdr <- p.adjust(wza1_means_outliers$emp_p,method = "fdr")
        wza2_means_outliers$fdr <- p.adjust(wza2_means_outliers$emp_p,method = "fdr")
        
        # # Visualise
        # ggplot(wza2_means,aes(x=blast_wza,y=outlier,fill=outlier))+
        #   stat_density_ridges(quantile_lines = TRUE)
        # 
        # ggplot(wza1_means_outliers,aes(x=emp_p))+
        #   geom_histogram(bins=100)+
        #   geom_hline(yintercept = (1-cutoff)*nrow(wza1_means_outliers),colour="red2")
        # 
        # ggplot(wza2_means_outliers,aes(x=emp_p))+
        #   geom_histogram(bins=100)+
        #   geom_hline(yintercept = (1-cutoff)*nrow(wza2_means_outliers),colour="red2")
        
        # # For each we want to identify the Higher Criticism Threshold...
        # wza1_hc <- HCthresh(wza1_means_outliers$emp_p,alpha=1,plotit = T)
        # wza2_hc <- HCthresh(wza2_means_outliers$emp_p,alpha=1-wza_cutoffs[1],plotit = T)
        
        # # Outlier detection based on Higher Ciritcism Threshold...
        # wza1_hc <- hc.thresh(wza1_means_outliers$emp_p,plot=T)
        # wza2_hc <- hc.thresh(wza2_means_outliers$emp_p,plot=T)
        # 
        # ggplot(wza1_means_outliers,aes(emp_p))+
        #   geom_histogram(binwidth=0.01)
        # 
        # # Collect the outliers...
        # wza1_final <- wza1_means_outliers[wza1_means_outliers$emp_p < wza1_hc,]
        # wza2_final <- wza2_means_outliers[wza2_means_outliers$emp_p < wza2_hc,]
        
        # Save the proportion that are true outliers...
        converge_p_matrix[row,col] <- length(wza1_means_outliers$fdr[wza1_means_outliers$fdr < 0.1])/nrow(wza1_means_outliers)
        converge_p_matrix[col,row] <- length(wza1_means_outliers$fdr[wza2_means_outliers$fdr < 0.1])/nrow(wza2_means_outliers)
        
        # },mc.cores=6)
        
        # # Extract results...
        # for(i in 1:nrow(converge_corr_matrix)){
        #   for(j in 1:ncol(converge_corr_matrix)){
        #     if(i != j){
        #       converge_corr_matrix[i,j] <- unlist(row_res[[i]][j])
        #     }
        #   }
        # }
        
      }
    } 
    # # Save these
    # saveRDS(converge_enrich_matrix_list,paste0("outputs/tmp_blast-aware_convergence_enrichment_matrices_",focal_climate,".rds"))
  }
  # Save these
  saveRDS(converge_p_matrix,paste0("outputs/tmp_blast-aware_convergence_p_matrix_",focal_climate,".rds"))
  
} else {
  res_mats <- readRDS(paste0("outputs/tmp_blast-aware_convergence_enrichment_matrices_",focal_climate,".rds"))
}
# converge_corr_matrix <- res_mats[[1]]
# converge_enrich_matrix <- res_mats[[2]]

# Add in the colnames rownames
for(i in 1:length(res_mats)){
  colnames(res_mats[[i]]) <- rownames(res_mats[[i]]) <- focal_datasets
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
    theme(axis.text.x = element_text(angle=45,hjust=1))+
    theme(axis.title = element_blank())
}

# Get OF tree
OF_tree <- ape::read.tree("outputs/orthology/Results_210913_17_genomes_with_conifers_noAA_filter/SpeciesTree_rooted.txt")

plot_matrix(res_mats[[2]])

test <- na.omit(reshape2::melt(res_mats[[2]]))
head(test[order(-test$value),])


# Compare to climate ------------------------------------------------------

# Fetch the climate data for all datasets...
all_climates <- data.frame(rbindlist(lapply(focal_datasets,function(dataset){
  climate_data <- read.table(paste0("outputs/GEA_res/",dataset,"/climate_cline.tsv"),header=T)
  climate_data$dataset <- dataset
  return(climate_data[,c("Long","Lat",focal_climate,"dataset")])
})))
colnames(all_climates)[3] <- "climate"

# Calculate NBR
NBR_res <- data.frame(dataset=focal_datasets,
                      NBR=NA)
for(dataset in NBR_res$dataset){
  NBR_res[NBR_res$dataset == dataset,"NBR"] <- log2(var(all_climates[all_climates$dataset==dataset,"climate"],na.rm = T)/var(all_climates$climate,na.rm = T))
}

# Calculate overlap-D
D_overlap_mat <- matrix(ncol=length(focal_datasets),nrow=length(focal_datasets))
colnames(D_overlap_mat) <- rownames(D_overlap_mat) <- focal_datasets
for(i in 1:nrow(D_overlap_mat)){
  for(j in 1:ncol(D_overlap_mat)){
    if(i < j){
      tmp_i <- na.omit(all_climates[all_climates$dataset == focal_datasets[i],"climate"])
      tmp_j <- na.omit(all_climates[all_climates$dataset == focal_datasets[j],"climate"])
      tmp_i_range <- round(seq(min(tmp_i),max(tmp_i)),2)
      tmp_j_range <- round(seq(min(tmp_j),max(tmp_j)),2)
      
      tmp_intersection <- intersect(tmp_i_range,tmp_j_range)
      if(length(tmp_intersection)!=0){
        D_overlap_mat[i,j] <- (max(intersect(tmp_i_range,tmp_j_range))-min(intersect(tmp_i_range,tmp_j_range)))/(max(c(tmp_i_range,tmp_j_range))-min(c(tmp_i_range,tmp_j_range)))
      } else {
        D_overlap_mat[i,j] <- 0
      }
    }
  }
}
plot_matrix(D_overlap_mat)

# Correlate...
D_melt <- na.omit(reshape2::melt(D_overlap_mat))
D_melt$comparison <- apply(D_melt[,c("Var1","Var2")],1,paste,collapse="-")
colnames(D_melt)[3] <- "climate_D"

enrich_melt <- na.omit(reshape2::melt(res_mats[[2]]))
enrich_melt$comparison <- apply(enrich_melt[,c("Var1","Var2")],1,paste,collapse="-")
colnames(enrich_melt)[3] <- "gene_enrich"

# Merge
enrich_melt <- merge(enrich_melt,D_melt[,c("climate_D","comparison")],by="comparison")
ggplot(enrich_melt,aes(climate_D,gene_enrich))+
  geom_point()


# Fetch branch lengths from genome comparisons...
library(adephylo)
tree_dists <- reshape2::melt(as.matrix(distTips(OF_tree)))
dataset_phylo_dists <- matrix(ncol=length(focal_datasets),nrow=length(focal_datasets))
rownames(dataset_phylo_dists) <- colnames(dataset_phylo_dists) <- focal_datasets
for(i in 1:nrow(dataset_phylo_dists)){
  for(j in 1:ncol(dataset_phylo_dists)){
    if(i < j ){
      genome_i <- as.character(na.omit(dataset_meta[dataset_meta$dataset == focal_datasets[i],"genome"]))
      genome_j <- as.character(na.omit(dataset_meta[dataset_meta$dataset == focal_datasets[j],"genome"]))
      
      dataset_phylo_dists[i,j] <- na.omit(c(as.matrix(distTips(OF_tree))[genome_i,genome_j],
                                            as.matrix(distTips(OF_tree))[genome_j,genome_i]))[1]
    }
  }
}
dataset_phylo_melt <- na.omit(reshape2::melt(dataset_phylo_dists))
dataset_phylo_melt$comparison <- apply(dataset_phylo_melt[,c("Var1","Var2")],1,paste,collapse="-")
colnames(dataset_phylo_melt)[3] <- "genome_distance"

# Merge and model
enrich_melt <- merge(enrich_melt,dataset_phylo_melt[,c("genome_distance","comparison")],by="comparison")
ggplot(enrich_melt,aes(genome_distance,gene_enrich))+
  geom_point()

# Try and build a crappy model
glm1 <- glm(log2(gene_enrich)~genome_distance*climate_D,data=enrich_melt,family="gaussian")
drop1(glm1,test="F")
glm2 <- glm(log2(gene_enrich)~genome_distance+climate_D,data=enrich_melt,family="gaussian")
drop1(glm2,test="F")

library(effects)
plot(allEffects(glm2))

# Visualise
enrich_melt$genome_groups <- cut_interval(enrich_melt$genome_distance,length = 0.2)
enrich_melt$climate_groups <- cut_interval(enrich_melt$climate_D,n = 4)

ggplot(enrich_melt,aes(genome_distance,log2(gene_enrich)))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~climate_groups,nrow = 1)

# mantels
library(vegan)
enrich_mat <- res_mats[[2]]
enrich_mat[lower.tri(enrich_mat)] <- enrich_mat[upper.tri(enrich_mat)]
D_overlap_mat[lower.tri(D_overlap_mat)] <- D_overlap_mat[upper.tri(D_overlap_mat)]
dataset_phylo_dists[lower.tri(dataset_phylo_dists)] <- dataset_phylo_dists[upper.tri(dataset_phylo_dists)]

mantel.partial(log2(enrich_mat), D_overlap_mat, dataset_phylo_dists, 
               method = "spearman", permutations = 999, na.rm = T, parallel = 6)

# # Fetch the ortho tree and plot them according to some vague phylo order...
# OF_tree <- ape::read.tree("outputs/orthology/Results_210913_17_genomes_with_conifers_noAA_filter/SpeciesTree_rooted.txt")
# phylo_order <- na.omit(dataset_meta)
# phylo_order <- data.frame(rbindlist(lapply(OF_tree$tip.label,function(x) phylo_order[phylo_order$genome == x,])))
# 
# test <- converge_enrich_matrix
# test[] <- 0
# heatmap(test)
# 
# # Look at these
# test <- reshape2::melt(converge_enrich_matrix)
# test <- test[order(-test$value),]
# head(test)
# test[test$value > 3,]
# table(test[test$value > 3,"Var1"])
# test[test$Var1 == "Amaranthus_tuberculatus_Wright_Individual",]

# # Attempt at Z-score based approach ---------------------------------------
# library(pbmcapply)
# Z_perms = 1000
# 
# wza1_outliers <- wza1[wza1$wza > quantile(wza1$wza,0.99),]
# wza2_outliers <- wza2[wza2$wza > quantile(wza2$wza,0.99),]
# 
# # It might be quicker to just expand out wza1_outliers once rather than repeatedly sample within the loops...
# blast_res_expanded <- data.frame(rbindlist(lapply(unique(wza1_outliers$qseqid),function(qseqid){
#   tmp <- blast_res_clean[blast_res_clean$qseqid==qseqid,]
#   tmp[rep(seq_len(nrow(tmp)), times = tmp$bitscore), ]
# })))
# 
# perm_res <- pbmclapply(1:Z_perms,function(iter){
#   set.seed(iter)
#   # # Based on bitscores, randomly draw 1 gene from species 2 based on the outliers from species 1
#   # tmp_ortho_fetch <- sapply(wza1_outliers$qseqid,function(qseqid){
#   #   sample(blast_res_clean[blast_res_clean$qseqid==qseqid,"sseqid"],prob = blast_res_clean[blast_res_clean$qseqid==qseqid,"bitscore"],1)
#   # })
#   
#   # Take a sample from the expanded blast_res_clean
#   sample_tmp <- blast_res_expanded[sample(1:nrow(blast_res_expanded)),]
#   sample_tmp <- sample_tmp[!(duplicated(sample_tmp$qseqid)),]
#   tmp_ortho_fetch <- sample_tmp$sseqid
#   
#   # How many of these are similarly in the top 1% - We expect 1%...
#   length(tmp_ortho_fetch[tmp_ortho_fetch %in% wza2_outliers$sseqid])/nrow(wza2_outliers)
# },mc.cores=6)
# 
# # Our Z score is then something along the lines of the median of this distribution compared with the null or something else?
# hist(unlist(perm_res))
# median(unlist(perm_res))/0.01

