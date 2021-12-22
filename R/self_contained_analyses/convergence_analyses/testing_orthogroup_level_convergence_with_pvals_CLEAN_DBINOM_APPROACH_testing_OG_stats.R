# Orthogroup-level association with environment...
lib <- c("mvmeta","qvalue","tidyr","ggtree","ape","VGAM","ggExtra","WeMix","pbmcapply","parallel","wCorr","data.table","ggplot2","viridis","ggridges","dplyr","patchwork","readr")
sapply(lib,library,character.only=T)

# Function Library --------------------------------------------------------

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
  return(tmp[tmp$gene_name != paste0(tmp$genome[1],"_"),])
},mc.cores = n_cores)))

# Merge it with OG_sequences
OG_sequences2 <- merge(OG_sequences,all_OG_long[,c("Orthogroup","gene_name")],by="gene_name")
OG_sequences2$full_OF <- OG_sequences2$OF_gene


# Group OGs into unique,  partial, full -----------------------------------
OG_genome <- unique(all_OG_long[,c("Orthogroup","genome")])
OG_counts <- data.frame(table(OG_genome$Orthogroup))
colnames(OG_counts) <- c("Orthogroup","species_count")
hist(OG_counts$species_count)

# Separate out partial and full based on species count of 10+
full_count=17
OG_counts[OG_counts$species_count == 1,"OG_type"] <- "Unique"
OG_counts[OG_counts$species_count > 1 &
            OG_counts$species_count < full_count,"OG_type"] <- "Partial"
OG_counts[OG_counts$species_count >= full_count,"OG_type"] <- "Full"

table(OG_counts$OG_type)

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
  
  # Replace emp pvals with sampled noise...
  # wza_merge$pvalue <- empPvals(stat=sample(wza_merge$weiZ),stat0=wza_merge$weiZ)
  
  # Mark
  wza_merge$dataset <- dataset
  wza_merge$genome <- genome
  
  # Add OG type to each data...
  wza_merge <- merge(wza_merge,OG_counts[,c("Orthogroup","OG_type")])
  
  return(data.table(wza_merge[,c("dataset","Orthogroup","pvalue","OG_type")]))
},mc.cores=6)
names(all_gea_res) <- focal_datasets


# Where do WZA outliers fall in each OG type ------------------------------
OG_type_wza_outliers <- rbindlist(all_gea_res)[,.(type_sum=nrow(.SD),type_signif=length(.SD$pvalue[.SD$pvalue < 0.05])),by=.(dataset,OG_type)]
gene_sums <- OG_type_wza_outliers[,.(gene_sum=sum(.SD$type_sum)),by=dataset]
OG_type_wza_outliers <- merge(OG_type_wza_outliers,gene_sums,by="dataset")

# Calculate expecteds...
OG_type_wza_outliers$obs_signif <- OG_type_wza_outliers$type_signif/OG_type_wza_outliers$type_sum
OG_type_wza_outliers$exp_signif <- (OG_type_wza_outliers$type_sum/OG_type_wza_outliers$gene_sum)*0.05*OG_type_wza_outliers$gene_sum

# Plot
to_plot <- melt(OG_type_wza_outliers[,.(dataset,OG_type,type_signif,exp_signif)])
ggplot(to_plot,aes(y=dataset,x=value,fill=variable))+
  geom_bar(stat="identity",position="dodge")+
  facet_wrap(~OG_type,ncol=1)

# Summary statistics for OGs ----------------------------------------------
# Quantify paralogy and species coverage...
paralog_counts <- rbindlist(all_gea_res)[,.(paralog_N=nrow(.SD)),by=.(Orthogroup,dataset)]
OG_stats <- paralog_counts[,.(mean_paralog_N=mean(.SD$paralog_N),
                              species_coverage=nrow(.SD)),by=Orthogroup]
hist(OG_stats$species_coverage)



# Orthogroup N significant dbinom pvals approach ---------------------------------
OG_maxP <- rbindlist(all_gea_res)[,.(maxP=min(.SD$pvalue),Ngenes_per_species=length(.SD$pvalue)),by=.(Orthogroup,dataset)]

# Choose a random set of OGs to work with...
OG_maxP <- OG_maxP[OG_maxP$Orthogroup %in% OG_maxP$Orthogroup[duplicated(OG_maxP$Orthogroup)],]
OG_subs <- OG_maxP$Orthogroup

OG_maxP <- OG_maxP[OG_maxP$Orthogroup %in% OG_subs,]

# Also just filter all the gea res
all_gea_res_sub <- rbindlist(lapply(all_gea_res,function(x) x[x$Orthogroup %in% OG_subs,]))

OG_maxP$Orthogroup_F <- factor(OG_maxP$Orthogroup,levels=names(table(OG_maxP$Orthogroup)))

# Get counts of OGs...
OG_counts <- table(OG_maxP$Orthogroup)

# Define cutoffs
cutoffs <- seq(0.001,0.05,0.001)

# Bonferroni
OG_maxP$maxP_bonf <- OG_maxP$maxP*OG_maxP$Ngenes_per_species
OG_maxP$maxP_bonf[OG_maxP$maxP_bonf > 1] <- 1

# Fetch cutoff-specific distributions of pvalues
pval_dists <- pbmclapply(cutoffs,function(flex_cutoff){
  # flex_cutoff <- 0.05
  
  # Define significant orthogroups based on multiple-testing corrections
  OG_maxP$maxP_bonf_signif  <- OG_maxP$maxP_bonf <= flex_cutoff
  
  # Sum up within all OGs the number of significant pvals
  OG_maxP_signif <- data.table(OG_maxP)[,.(NSignif=sum(.SD$maxP_bonf_signif),N=length(.SD$maxP_bonf_signif)),by=Orthogroup]
  
  # We can now just remove any cases that we don't want to test where NSignif < 2
  OG_maxP_signif <- OG_maxP_signif[OG_maxP_signif$NSignif > 1,]
  OG_maxP_signif$final_pval <- sapply(1:nrow(OG_maxP_signif),function(x) {
    sum(dbinom(OG_maxP_signif$NSignif[x]:OG_maxP_signif$N[x],OG_maxP_signif$N[x],prob = flex_cutoff))
  })
  
  out <- data.frame(Orthogroup=OG_maxP_signif$Orthogroup,
                    pvalue=OG_maxP_signif$final_pval,
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
final_score

# FINAL OG SCORES...
final_OG_res <- data.frame(rbindlist(lapply(pval_dists,'[[',1)))
final_OG_res_scores <- sort(table(final_OG_res[final_OG_res$fdr < 0.05,"Orthogroup"]))


# Compare results against orthogroup statistics ---------------------------
final_OG_res_scores_merge <- merge(final_OG_res[final_OG_res$cutoff==0.05,],OG_stats,by="Orthogroup")
final_OG_res_scores_merge$pvalue_group <- cut_interval(final_OG_res_scores_merge$pvalue,n = 20)
final_OG_res_scores_merge$fdr_group <- cut_interval(final_OG_res_scores_merge$fdr,n = 20)

# Compare with rate of paralogy...
ggplot(final_OG_res_scores_merge,aes(x=species_coverage,y=fdr_group))+
  geom_boxplot()

# Explore cases where we saw a big shift between low to high pval
OG_maxP$bonf_diff <- OG_maxP$maxP_bonf - OG_maxP$maxP
head(OG_maxP[order(-OG_maxP$bonf_diff),])
