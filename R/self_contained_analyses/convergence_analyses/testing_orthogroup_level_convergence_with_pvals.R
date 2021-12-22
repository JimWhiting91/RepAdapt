# Orthogroup-level association with environment...
lib <- c("mvmeta","qvalue","tidyr","ggtree","ape","VGAM","ggExtra","WeMix","pbmcapply","parallel","wCorr","data.table","ggplot2","viridis","ggridges","dplyr","patchwork","readr")
sapply(lib,library,character.only=T)


# Function Library --------------------------------------------------------
pMaxTest <- function( p, setSize){
  numTests = length(p)
  sortedP = sort(p, decreasing = T)
  sapply(2:setSize, function(setSize_x) {
    pThreshold = sortedP[1 + numTests - setSize_x] 
    sum( dbinom( setSize_x:numTests, numTests, pThreshold))
  })
}

# permute_OG_stats <- function(full_OG_data,focal_OG,nperms,pval_column="pval_column"){
#   
#   perm_res <- lapply(1:nperms,function(perm){
#     
#     # Randomise
#     tmp <- full_OG_data
#     tmp$Orthogroup <- sample(tmp$Orthogroup,replace = F)
#     tmp <- tmp[tmp$Orthogroup == focal_OG,]
#     
#     # Return per permutation stats 
#     return(list(sum=sum(-log10(tmp[,pval_column])),
#                 Nsignif=length(tmp[,pval_column][tmp[,pval_column] < 0.05]),
#                 sum_NoTop=sum(-log10(tmp[,pval_column][tmp[,pval_column] != min(tmp[,pval_column])]))))
#   })
#   
#   out <- data.frame(sum=sapply(perm_res,'[[',1),
#                     Nsignif=sapply(perm_res,'[[',2),
#                     sum_NoTop=sapply(perm_res,'[[',1))
#   return(out)
# }

# Calculate covergence stats per OG
calculate_OG_stats <- function(OG_data,
                               pval_column="maxP",
                               dataset_id_column="dataset"){
  
  pval_vector <- data.frame(OG_data[,pval_column])[,1]
  names(pval_vector) <- data.frame(OG_data[,dataset_id_column])[,1]
  
  # Return per permutation stats 
  return(list(summaries=data.frame(sum=sum(-log10(OG_data[,pval_column])),
                                   Nsignif=length(OG_data[,pval_column][OG_data[,pval_column] < 0.05]),
                                   sum_NoTop=sum(-log10(OG_data[,pval_column][OG_data[,pval_column] != min(OG_data[,pval_column])]))),
              pvals=pval_vector))
  
}
# Calculate covergence stats per OG
# calculate_OG_stats_DT <- function(OG_data,
#                                pval_column="maxP",
#                                dataset_id_column="dataset"){
#   
#   pval_vector <- OG_data[,get(pval_column)]
#   names(pval_vector) <- OG_data[,get(dataset_id_column)]
#   
#   # # Return per permutation stats 
#   # return(list(summaries=data.frame(sum=sum(-log10(pval_vector)),
#   #                                  Nsignif=length(pval_vector[pval_vector < 0.05]),
#   #                                  sum_NoTop=sum(-log10(pval_vector[pval_vector != min(pval_vector)]))),
#   #             pvals=pval_vector))
#   # Return per permutation stats
#   return(list(summaries=data.frame(sum=sum(-log10(pval_vector)),
#                                    Nsignif=length(pval_vector[pval_vector < 0.05]),
#                                    sum_NoTop=sum(-log10(pval_vector[pval_vector != min(pval_vector)]))),
#               pvals=pval_vector))
#   
# }


calculate_OG_stats_DT <- function(pval_vector){
  
  return(summaries=data.frame(sum=sum(-log10(pval_vector)),
                              Nsignif=length(pval_vector[pval_vector < 0.05]),
                              sum_NoTop=sum(-log10(pval_vector[pval_vector != min(pval_vector)]))))
}



# calculate_OG_stats_plus <- function(OG_data,
#                                pval_column="maxP",
#                                dataset_id_column="dataset",
#                                expected_pvals){
#   
#   # Take pvalues like normal
#   pval_vector <- data.frame(OG_data[,pval_column])[,1]
#   names(pval_vector) <- data.frame(OG_data[,dataset_id_column])[,1]
#   
#   # Correct them based on the expecteds drawn from all permutations
#   pval_vector_corrected <- sapply(1:length(pval_vector),function(x){
#     empPvals(-log10(pval_vector[x]),stat0 = -log10(expected_pvals[x,]))
#   })
#   
#   # bpval_vector <- data.frame(OG_data[,"bonfP"])[,1]
#   # names(bpval_vector) <- data.frame(OG_data[,dataset_id_column])[,1]
#   
#   # Return per permutation stats 
#   return(list(summaries=data.frame(sum=sum(-log10(OG_data[,pval_column])),
#                                    Nsignif=length(OG_data[,pval_column][OG_data[,pval_column] < 0.05]),
#                                    sum_NoTop=sum(-log10(OG_data[,pval_column][OG_data[,pval_column] != min(OG_data[,pval_column])])),
#                                    Nsignif_plus=length(pval_vector_corrected[pval_vector_corrected < 0.05])),
#               pvals=pval_vector))
#   # bonf=bpval_vector))
# }
############################################################################
focal_climate <- "mean_diurnal"

n_cores=4
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
  
  # Calclate emp pvals
  wza_merge$pvalue <- empPvals(stat=wza_merge$weiZ,stat0=wza_merge$weiZ)
  # wza_merge$pvalue <- pnorm(wza_merge$weiZ, lower.tail=FALSE)
  # wza_merge$pvalue <- runif(nrow(wza_merge),0,1)
  
  # Mark
  wza_merge$dataset <- dataset
  wza_merge$genome <- genome
  return(data.table(wza_merge[,c("dataset","Orthogroup","pvalue")]))
},mc.cores=6)
names(all_gea_res) <- focal_datasets

# # Orthogroup Pmax approach ----------------------------------------
# OG_maxP <- data.frame(data.frame(rbindlist(all_gea_res)) %>%
#                         group_by(dataset,Orthogroup) %>%
#                         summarise(maxP=min(pvalue)))
# OG_maxP <- OG_maxP[OG_maxP$Orthogroup %in% OG_maxP$Orthogroup[duplicated(OG_maxP$Orthogroup)],]
# OG_maxP$Orthogroup_F <- factor(OG_maxP$Orthogroup,levels=names(table(OG_maxP$Orthogroup)))
# 
# OG_pMax_vectors <- OG_maxP %>% group_by(Orthogroup_F) %>%
#   group_map(~pMaxTest(.x$maxP,length(.x$maxP)))
# names(OG_pMax_vectors) <- levels(OG_maxP$Orthogroup_F)
# 
# # Sanity check
# table(sapply(OG_pMax_vectors,length) == table(OG_maxP$Orthogroup)-1)
# 
# # Build a reference of idnex to OG name
# OG_index <- data.frame(OG=names(OG_pMax_vectors),
#                        index=1:length(names(OG_pMax_vectors) ))
# 
# OG_pMax_dd <- data.frame(pvals=unlist(OG_pMax_vectors))
# OG_pMax_dd$OG_index <- rep(paste0("OG_",1:length(OG_pMax_vectors)),times=sapply(OG_pMax_vectors,length))
# OG_lengths <- table(OG_pMax_dd$OG_index)+1
# OG_lengths <- OG_lengths[unique(OG_pMax_dd$OG_index)]
# OG_pMax_dd$setSize <- unlist(lapply(OG_lengths,function(i) seq(2,i,1)))
# 
# # Merge set size and OG
# OG_pMax_dd$merge_var <- paste0(OG_pMax_dd$OG_index,":",OG_pMax_dd$setSize)
# 
# permN=100
# perm_res <- pbmclapply(1:permN,function(perm){
#   set.seed(perm)
# 
#   # First shuffle
#   all_gea_res_shuff <- lapply(all_gea_res,function(tmp){
#     tmp$Orthogroup <- sample(tmp$Orthogroup)
#     return(tmp)
#   })
# 
#   # Now do pmaxs
#   OG_maxP_shuf <- data.frame(data.frame(rbindlist(all_gea_res_shuff)) %>%
#                                group_by(dataset,Orthogroup) %>%
#                                summarise(maxP=min(pvalue)))
# 
#   OG_maxP_shuf <- OG_maxP_shuf[OG_maxP_shuf$Orthogroup %in% OG_maxP_shuf$Orthogroup[duplicated(OG_maxP_shuf$Orthogroup)],]
#   OG_maxP_shuf$Orthogroup_F <- factor(OG_maxP_shuf$Orthogroup,levels=names(table(OG_maxP_shuf$Orthogroup)))
# 
#   OG_pMax_vectors_shuf <- OG_maxP_shuf %>% group_by(Orthogroup_F) %>%
#     group_map(~pMaxTest(.x$maxP,length(.x$maxP)))
# 
#   # Reformat the data...
#   out <- data.frame(pvals=unlist(OG_pMax_vectors_shuf))
#   out$OG_index <- rep(paste0("OG_",1:length(OG_pMax_vectors)),times=sapply(OG_pMax_vectors,length))
#   OG_lengths <- table(out$OG_index)+1
#   OG_lengths <- OG_lengths[unique(out$OG_index)]
#   out$setSize <- unlist(lapply(OG_lengths,function(i) seq(2,i,1)))
# 
#   return(out)
# },mc.cores=n_cores)
# 
# # Fetch all cutoffs
# permuted_OG_cutoffs <- data.frame(rbindlist(perm_res) %>% group_by(OG_index,setSize) %>%
#   summarise(cutoff=quantile(pvals,probs=0.05)))
# 
# # Merge set size and OG
# permuted_OG_cutoffs$merge_var <- paste0(permuted_OG_cutoffs$OG_index,":",permuted_OG_cutoffs$setSize)
# 
# # Merge observed with cutoffs...
# OG_pMax_merge <- merge(OG_pMax_dd,permuted_OG_cutoffs[,c("cutoff","merge_var")])
# 
# nrow(OG_pMax_merge[OG_pMax_merge$pvals < OG_pMax_merge$cutoff,])/nrow(OG_pMax_merge)
# 
# # Retain 'significant' OGs...
# signif_pmax <- OG_pMax_merge[OG_pMax_merge$pvals < OG_pMax_merge$cutoff,]
# 
# # How many sets are significant per OG
# signif_OG_counts <- table(signif_pmax$OG_index)
# 
# # How many OGs are significant at all...
# table(unique(OG_pMax_merge$OG_index) %in% unique(signif_pmax$OG_index))
# 
# OG_maxP[OG_maxP$Orthogroup==OG_index$OG[7333],]
# 
# permued_OG_cutoffs[[2]]


# Orthogroup N significant pvals approach ---------------------------------
# OG_maxP <- data.frame(data.frame(rbindlist(all_gea_res)) %>%
#                         group_by(dataset,Orthogroup) %>%
#                         summarise(wza=max(weiZ),
#                                   maxP=min(pvalue),
#                                   bonfP=maxP*length(pvalue)))
# OG_maxP <- data.frame(rbindlist(all_gea_res) %>%
#                         group_by(dataset,Orthogroup) %>%
#                         summarise(maxP=min(pvalue)))

OG_maxP <- rbindlist(all_gea_res)[,.(min(.SD$pvalue),length(.SD$pvalue)),by=.(Orthogroup,dataset)]
colnames(OG_maxP) <- c("Orthogroup","dataset","maxP","Ngenes_per_species")

# Choose a random set of OGs to work with...
OG_maxP <- OG_maxP[OG_maxP$Orthogroup %in% OG_maxP$Orthogroup[duplicated(OG_maxP$Orthogroup)],]
# OG_subs <- sort(sample(OG_maxP$Orthogroup,1000))
OG_subs <- OG_maxP$Orthogroup
# OG_subs <- rbindlist(all_gea_res)[,.(head(sort(.SD$pvalue),100),.SD$Orthogroup),by=.(dataset)]

OG_maxP <- OG_maxP[OG_maxP$Orthogroup %in% OG_subs,]

# Also just filter all the gea res
all_gea_res_sub <- rbindlist(lapply(all_gea_res,function(x) x[x$Orthogroup %in% OG_subs,]))

OG_maxP$Orthogroup_F <- factor(OG_maxP$Orthogroup,levels=names(table(OG_maxP$Orthogroup)))

# OG_index <- data.frame(Orthogroup=levels(OG_maxP$Orthogroup_F),
#                        index=1:length(levels(OG_maxP$Orthogroup_F)))

# Get counts of OGs...
OG_counts <- table(OG_maxP$Orthogroup)

# # Calculate observed stats
# OG_obvs_res <- OG_maxP %>% group_by(Orthogroup_F) %>%
#   group_map(~calculate_OG_stats(.x,"maxP"))
# names(OG_obvs_res) <- levels(OG_maxP$Orthogroup_F)

# Permute over all orthogroups, and pull number of low pvals...
permN=100

# # Shuffle all of our datasets...
# gea_res_shuff_perms <- pbmclapply(1:permN,function(perm){
# 
#   return(gea_res_shuff)
# },mc.cores=n_cores)


# ########################################################################################
# # Only store the permuted OG vectors
# shuffled_pvals <- setDT(pbmclapply(1:permN,function(perm){
#   set.seed(perm)
#   
#   # Resample the pvals
#   tmp <- all_gea_res_sub[,.(sample(.SD$pvalue,replace = F)),by="dataset"]
#   tmp$Orthogroup <- all_gea_res_sub$Orthogroup
#   
#   # Return a matrix of all the resampled pvalues
#   tmp[,min(.SD$V1),by=.(Orthogroup,dataset)]$V1
# },mc.cores=6))
# ########################################################################################

# Define cutoffs
cutoffs <- seq(0.001,0.05,0.001)

# Bonferroni correct
# OG_maxP$maxP_bonf <- OG_maxP$maxP*OG_maxP$Ngenes_per_species
# Dunn-Sidak correction
OG_maxP$maxP_DS <- 1 - (1 - OG_maxP$maxP)^OG_maxP$Ngenes_per_species

# OG_pMax_vectors <- OG_maxP %>% group_by(Orthogroup) %>%
#   group_map(~pMaxTest(.x$maxP_DS,length(.x$maxP_DS)))

# Fetch cutoff-specific distributions of pvalues
pval_dists <- pbmclapply(cutoffs,function(flex_cutoff){
  # flex_cutoff <- 0.05
  
  # Define significant orthogroups based on multiple-testing corrections
  OG_maxP$maxP_DS_signif  <- OG_maxP$maxP_DS <= flex_cutoff
  
  # # gene N FPR
  # fpr_res <- data.frame(data.table(OG_maxP)[,.(sum(.SD$maxP_signif)/nrow(.SD)),by="Ngenes_per_species"])
  # ggplot(fpr_res,aes(Ngenes_per_species,V1))+
  #   geom_point()
  
  # Sum up within all OGs the number of significant pvals
  OG_maxP_signif <- data.table(OG_maxP)[,.(sum(.SD$maxP_DS_signif),length(.SD$maxP_DS_signif)),by=Orthogroup]
  colnames(OG_maxP_signif) <- c("Orthogroup","NSignif","N")
  OG_maxP_signif$final_pval <- sapply(1:nrow(OG_maxP_signif),function(x) sum(dbinom(OG_maxP_signif$NSignif[x]:OG_maxP_signif$N[x],OG_maxP_signif$N[x],prob = flex_cutoff)))
  
  out <- data.frame(Orthogroup=OG_maxP_signif$Orthogroup,
                    fdr=p.adjust(OG_maxP_signif$final_pval,method="fdr"),
                    cutoff=flex_cutoff)
  return(list(orthogroup_res=out,
              nsignif=length(out$fdr[out$fdr < 0.05])))
},mc.cores=n_cores)

to_plot <- data.frame(cutoffs=cutoffs,
                      signif_OG=sapply(pval_dists,'[[',2))
ggplot(to_plot,aes(cutoffs,signif_OG))+
  geom_point()

# FINAL DATASET SCORE
final_score <- mean(to_plot$signif_OG)

# FINAL OG SCORES...
final_OG_res <- data.frame(rbindlist(lapply(pval_dists,'[[',1)))
final_OG_res_scores <- sort(table(final_OG_res[final_OG_res$fdr < 0.05,"Orthogroup"]))

# How many significant orthogroups...
length(final_OG_res_scores)

# pMAX version of test ----------------------------------------------------
# Calculate pMaxTest for all genes...
OG_pMax_res <- OG_maxP[,.(min(pMaxTest(.SD$maxP_DS,length(.SD$maxP_DS))),
                          which(pMaxTest(.SD$maxP_DS,length(.SD$maxP_DS))==min(pMaxTest(.SD$maxP_DS,length(.SD$maxP_DS)))),
                          length(pMaxTest(.SD$maxP_DS,length(.SD$maxP_DS)))),
                       by=Orthogroup]
colnames(OG_pMax_res) <- c("Orthogroup","pMax_p","NSignif","NumTests")

# Multiple correct again...
OG_pMax_res$pMax_p_DS <- 1 - (1 - OG_pMax_res$pMax_p)^OG_pMax_res$NumTests
OG_pMax_res$fdr <- p.adjust(OG_pMax_res$pMax_p_DS,method = "fdr")
final_pmax_res <- OG_pMax_res[OG_pMax_res$fdr < 0.05,]

names(final_OG_res_scores) %in% final_pmax_res$Orthogroup
names(final_OG_res_scores)[!(names(final_OG_res_scores) %in% final_pmax_res$Orthogroup)]
# OG_obvs_res <- data.table(OG_maxP)[,calculate_OG_stats_DT(.SD$maxP),by=Orthogroup_F]
# all_pvals <- split(data.table(OG_maxP),by="Orthogroup_F")
# 
# 
# shuffled_pvals$Orthogroup <- all_gea_res_sub$Orthogroup
# shuffled_pvals$dataset <- all_gea_res_sub$dataset
# 
# # Run the function over all 
# test <- pbmclapply(1:permN,function(perm){
#   shuffled_pvals[,calculate_OG_stats_DT(.SD),
#                                               by=.(dataset,Orthogroup)]
# 
# # Permute slowly
# perm_res <- pbmclapply(1:permN,function(perm){
#   print(perm)
#   # set.seed(perm)
#   
#   # Set up new OGs
#   all_gea_res_sub$rand_OG <- shuffled_OGs[[perm]]
#   
#   # Re-estimate max Wza
#   OG_maxP_shuff <- all_gea_res_sub[,min(.SD$pvalue),by=.(rand_OG,dataset)]
#   colnames(OG_maxP_shuff)[ncol(OG_maxP_shuff)] <- "maxP"
#   # 
#   # OG_maxP_shuff <- data.frame(data.frame(rbindlist(gea_res_shuff)) %>%
#   #                         group_by(dataset,Orthogroup) %>%
#   #                         summarise(wza=max(weiZ),
#   #                                   maxP=min(pvalue),
#   #                                   bonfP=maxP*length(pvalue)))
#   
#   
#   # Fetch permutations stats for all OGs
#   # Data table equivalent...
#   OG_obvs_res <- OG_maxP_shuff[,calculate_OG_stats_DT(.SD$maxP),by=rand_OG]
#   all_pvals <- split(data.table(OG_maxP_shuff),by="rand_OG")
#   # OG_perm_res <- OG_maxP_shuff %>% group_by(Orthogroup_F) %>%
#   #   group_map(~calculate_OG_stats(.x,"maxP"))
#   # names(OG_perm_res) <- levels(OG_maxP$Orthogroup_F)
#   
#   return(list(OG_obvs_res,all_pvals))
# },mc.cores=6)
# 
# # Merge all of the OG permutations into a single df per OG
# perm_OG_summaries <- pbmclapply(1:length(OG_obvs_res),function(perm){
#   OG_perm_tmp <- lapply(perm_res,'[[',perm)
#   data.frame(rbindlist(lapply(OG_perm_tmp,'[[',1)))
# },mc.cores=n_cores)
# 
# # Build matrices of permuted maximum pvals
# perm_OG_maxP <- pbmclapply(1:length(OG_obvs_res),function(perm){
#   OG_perm_tmp <- lapply(perm_res,'[[',perm)
#   pval_matrix <- matrix(nrow=length(OG_perm_tmp[[1]][[2]]),unlist(lapply(OG_perm_tmp,'[[',2)))
#   rownames(pval_matrix) <- names(OG_perm_tmp[[1]][[2]])
#   pval_matrix
# },mc.cores=n_cores)
# names(perm_OG_maxP) <- names(OG_obvs_res)
# 
# # #### Perm ROUND 2 ####
# # perm_res2 <- pbmclapply(1:permN,function(perm){
# #   print(perm)
# #   set.seed(perm)
# #   
# #   # Shuffle within dataset
# #   gea_res_shuff <- lapply(all_gea_res,function(gea_res){
# #     gea_res <- gea_res[gea_res$Orthogroup %in% OG_subs,]
# #     gea_res$Orthogroup <- sample(gea_res$Orthogroup)
# #     gea_res
# #   })
# #   
# #   # Re-estimate max Wza
# #   OG_maxP_shuff <- data.frame(data.frame(rbindlist(gea_res_shuff)) %>%
# #                                 group_by(dataset,Orthogroup) %>%
# #                                 summarise(wza=max(weiZ),
# #                                           maxP=min(pvalue),
# #                                           bonfP=maxP*length(pvalue)))
# #   
# #   # Factorise
# #   OG_maxP_shuff$Orthogroup_F <- factor(OG_maxP_shuff$Orthogroup,levels=names(table(OG_maxP$Orthogroup)))
# #   
# #   # Fetch permutations stats for all OGs
# #   OG_perm_res <- OG_maxP_shuff %>% group_by(Orthogroup_F) %>%
# #     group_map(~calculate_OG_stats_plus(.x,pval_column = "maxP",expected_pvals = perm_OG_maxP[unique(.x$Orthogroup)]))
# #   names(OG_perm_res) <- levels(OG_maxP$Orthogroup_F)
# #   
# #   return(OG_perm_res)
# # },mc.cores=n_cores)
# # ######################
# # 
# # Now run over all observed values and calculate pvals for the 3 variables...
# OG_perm_pvals <- data.frame(rbindlist(pbmclapply(1:length(OG_obvs_res),function(OG){
#   
#   perm_p_res <- sapply(1:length(OG_obvs_res[[OG]]$pvals),function(species){
#     perm_p <- perm_OG_maxP[[OG]][species,]
#     empPvals(-log10(OG_obvs_res[[OG]]$pvals[species]),stat0 = -log10(perm_p))
#   })
#   Nsignif_plus_obs <- length(perm_p_res[perm_p_res <= 0.05])
#   
#   # # Transform our permuted pvalue matrix to empiricals
#   # empirical_perm_p <- perm_OG_maxP[[OG]]
#   # for(i in 1:nrow(empirical_perm_p)){
#   #   empirical_perm_p[i,] <- empPvals(-log10(perm_OG_maxP[[OG]][i,]),stat0=-log10(perm_OG_maxP[[OG]][i,]))
#   # }
#   # 
#   # # Go through empirical perm p matrix column-wise, and count <0.05
#   # col_perms <- sapply(1:ncol(empirical_perm_p),function(col){
#   #   length(empirical_perm_p[,col][empirical_perm_p[,col] < 0.05])
#   # })
#   
#   # Calculate pvals on own 95th cutoff
#   data.frame(sum_pval=empPvals(OG_obvs_res[[OG]]$summaries[,1],perm_OG_summaries[[OG]][,1]),
#              Nsignif_pval=empPvals(OG_obvs_res[[OG]]$summaries[,2],perm_OG_summaries[[OG]][,2]),
#              sum_NoTop_pval=empPvals(OG_obvs_res[[OG]]$summaries[,3],perm_OG_summaries[[OG]][,3]),
#              Nsignif_plus_pval=sum(dbinom(Nsignif_plus_obs:length(perm_p_res),length(perm_p_res),0.05)),
#              # Nsignif_plus_pval_stupid=empPvals(Nsignif_plus_obs,stat0=col_perms),
#              OG=OG_index$Orthogroup[OG])
#   
# },mc.cores = n_cores)))
# 
# # Plot them out...
# OG_perm_pvals$OG_count_n <- OG_counts
# hist(OG_perm_pvals$sum_pval)
# hist(OG_perm_pvals$Nsignif_pval)
# hist(OG_perm_pvals$sum_NoTop_pval)
# hist(OG_perm_pvals$Nsignif_plus_pval)
# 
# # Associations with other vars
# ggplot(OG_perm_pvals,aes(-log10(sum_pval),-log10(sum_NoTop_pval)))+
#   geom_point()+
#   stat_density_2d_filled(alpha=0.5)
# ggplot(OG_perm_pvals,aes(-log10(sum_pval),-log10(Nsignif_pval)))+
#   geom_point()+
#   stat_density_2d_filled(alpha=0.5)
# ggplot(OG_perm_pvals,aes(-log10(sum_NoTop_pval),-log10(Nsignif_pval)))+
#   geom_point()+
#   stat_density_2d_filled(alpha=0.5)
# 
# ggplot(OG_perm_pvals,aes(Nsignif_OnSelf_prop,-log10(Nsignif_pval)))+
#   geom_point()+
#   stat_density_2d_filled(alpha=0.5)
# 
# ggplot(OG_perm_pvals,aes(-log10(Nsignif_plus_pval),-log10(Nsignif_pval)))+
#   geom_point()+
#   stat_density_2d_filled(alpha=0.5)
# 
# # What is the effect of species N
# ggplot(OG_perm_pvals,aes(OG_count_n,-log10(sum_pval)))+
#   geom_point()+
#   stat_density_2d_filled(alpha=0.5)
# ggplot(OG_perm_pvals,aes(OG_count_n,-log10(sum_NoTop_pval)))+
#   geom_point()+
#   stat_density_2d_filled(alpha=0.5)
# ggplot(OG_perm_pvals,aes(OG_count_n,-log10(Nsignif_pval)))+
#   geom_point()+
#   stat_density_2d_filled(alpha=0.5)
# 
# # What are the top OGs...
# head(OG_perm_pvals[order(OG_perm_pvals$sum_pval),])
# head(OG_perm_pvals[order(OG_perm_pvals$sum_NoTop_pval),])
# head(OG_perm_pvals[order(OG_perm_pvals$Nsignif_pval),])
# head(OG_perm_pvals[order(OG_perm_pvals$Nsignif_plus_pval),])
# 
# nrow(OG_perm_pvals[OG_perm_pvals$sum_pval < 0.05,])/nrow(OG_perm_pvals)
# nrow(OG_perm_pvals[OG_perm_pvals$Nsignif_pval < 0.05,])/nrow(OG_perm_pvals)
# nrow(OG_perm_pvals[OG_perm_pvals$Nsignif_pval < 0.05 &
#                      OG_perm_pvals$sum_pval < 0.05,])/nrow(OG_perm_pvals)
# 
# # Check ranking effect...
# OG_perm_pvals$OG_index <- 1:nrow(OG_perm_pvals)
# ggplot(OG_perm_pvals,aes(OG_index,-log10(sum_pval)))+
#   geom_point()+
#   stat_density_2d_filled(alpha=0.5)
# hist(OG_perm_pvals[OG_perm_pvals$sum_pval < 0.05,"OG_index"])
# 




