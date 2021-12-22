# Compare pSONIC vs raw Orthofinder
lib <- c("ggridges","data.table","ggplot2","dplyr")
lapply(lib,library,character.only=T)


# Function library --------------------------------------------------------
count_singletons <- function(orthogroups){
  
  # Output mat
  out_mat <- matrix(nrow=ncol(orthogroups)-1,ncol=3)
  rownames(out_mat) <- colnames(orthogroups)[2:ncol(orthogroups)]
  
  # Populate
  for(spec in rownames(out_mat)){
    set_size <- sapply(orthogroups[,spec],function(x) length(strsplit(x,",")[[1]]))
    out_mat[spec,1] <- length(set_size[set_size==1])
    out_mat[spec,2] <- length(set_size[set_size==1])/length(set_size)
    out_mat[spec,3] <- sum(set_size)
  }
  return(out_mat)
}

count_number_genes_per_set <- function(orthogroups){
  
  # Output mat
  specs <- colnames(orthogroups)[colnames(orthogroups) != "Orthogroup"]
  out_mat <- matrix(nrow=nrow(orthogroups),ncol=length(specs))
  colnames(out_mat) <- specs
  
  # Populate
  for(spec in specs){
    out_mat[,spec] <- sapply(orthogroups[,spec],function(x) length(strsplit(x,",")[[1]]))
  }
  return(out_mat)
}

all_row_same <- function(row,same_value=1){
  all(row == same_value)
}

count_full_singleton_sets <- function(orthogroups){
  
  # Output mat
  specs <- colnames(orthogroups)[colnames(orthogroups) != "Orthogroup"]
  out_mat <- matrix(nrow=nrow(orthogroups),ncol=length(specs))
  colnames(out_mat) <- specs
  
  # Populate
  for(spec in specs){
    out_mat[,spec] <- sapply(orthogroups[,spec],function(x) length(strsplit(x,",")[[1]]))
  }
  
  # Count full singletons...
  singleton_rows <- apply(out_mat,1,all_row_same,same_value=1)
  
  return(data.frame(table(singleton_rows)))
}

# Analysis ----------------------------------------------------------------
# Where are the results?
results_dir <- "outputs/orthology/pSONIC_brassica/"

# Make results dir
dir.create("figs/orthofinder_res",showWarnings = F)
fig_dir <- "figs/orthofinder_res"

# Fetch both results
OF_groups <- data.frame(read.table(paste0(results_dir,"Orthogroups.tsv"),sep="\t",header=T))
pSONIC_groups <- data.frame(read.table(paste0(results_dir,"pSONIC.txt"),sep="\t",header=F))[,1:ncol(OF_groups)]
colnames(pSONIC_groups) <- colnames(OF_groups)


# Compare singletons ------------------------------------------------------
OF_sing <- count_singletons(OF_groups)
pSONIC_sing <- count_singletons(pSONIC_groups)

# Plot together
singleton_plot <- data.frame(rbind(OF_sing,pSONIC_sing))
colnames(singleton_plot) <- c("singleton_N","singleton_prop","total_genes_assigned")
singleton_plot$analysis <- rep(c("OF","pSONIC"),each=nrow(OF_sing))
singleton_plot$species <- rep(rownames(OF_sing),2)

# Singleton Count
sing_count <- ggplot(singleton_plot,aes(x=species,y=singleton_N,fill=analysis))+
  geom_bar(stat = "identity",position="dodge")+
  ggtitle("Comparison of Singletons")

# Proportion assigned singletons
prop_sing_assigned <- ggplot(singleton_plot,aes(x=species,y=singleton_N/total_genes_assigned,fill=analysis))+
  geom_bar(stat = "identity",position="dodge")+
  ggtitle("Proportion assigned singletons")

# Total genes assigned to OF
total_sing_assigned <- ggplot(singleton_plot,aes(x=species,y=total_genes_assigned,fill=analysis))+
  geom_bar(stat = "identity",position="dodge")+
  ggtitle("Total genes assigned to orthogroups")

# What is total number of orthogroups with 1 gene for all species...
pSONIC_full_singletons <- count_full_singleton_sets(pSONIC_groups)
OF_full_singletons <-  count_full_singleton_sets(OF_groups)
plot_full_singletons <- data.frame(props=c(pSONIC_full_singletons$Freq[2]/sum(pSONIC_full_singletons$Freq),
                                           OF_full_singletons$Freq[2]/sum(OF_full_singletons$Freq)),
                                   analysis=c("pSONIC","OF"))
full_singletons <- ggplot(plot_full_singletons,aes(x=analysis,y=props))+
  geom_bar(stat="identity")+
  ggtitle("Proportion of OG full singletons")

library(patchwork)
(sing_count/prop_sing_assigned)|(total_sing_assigned/full_singletons)

# Count genes in orthogroups in each --------------------------------------
pSONIC_counts <- count_number_genes_per_set(pSONIC_groups)
OF_counts <- count_number_genes_per_set(OF_groups)
pSONIC_count_melt <- reshape2::melt(pSONIC_counts)
pSONIC_count_melt$analysis <- "pSONIC"
OF_count_melt <- reshape2::melt(OF_counts)
OF_count_melt$analysis <- "OF"
to_plot <- rbind(pSONIC_count_melt,OF_count_melt)

ggplot(to_plot[to_plot$value < 20,],aes(x=value))+
  geom_histogram()+
  facet_grid(Var2~analysis)+
  #geom_bar(stat="identity")+
  ggtitle("Number of genes in orthogroup")


# Assess dropout from raw collinearity groups -----------------------------
raw_collinearity_groups <- read.table(paste0(results_dir,"pSONIC.RawGroups.txt"),sep="\t")
colnames(raw_collinearity_groups) <- c("group_length","PASS_genepairs","NOTPASS_genepairs","NoCall_genepairs","pattern","orientation")

# Get collinearity names...
raw_blocks <- read.table(paste0(results_dir,"pSONIC_brassica.collinearity"),comment.char = "#",fill=T)
colnames(raw_blocks) <- c("block","gene_pair","geneA","geneB","score")

# Add some dataset identifier to the groups summary...
raw_collinearity_groups$species_pair <- NA
for(i in 1:nrow(raw_collinearity_groups)){
  print(i)
  tmp <- raw_blocks[raw_blocks$block == paste0(i-1,"-"),]
  speciesA <- strsplit(tmp$geneA[1],"_")[[1]][1]
  speciesB <- strsplit(tmp$geneB[1],"_")[[1]][1]
  both_species <- sort(c(speciesA,speciesB))
  raw_collinearity_groups$species_pair[i] <- paste(both_species,collapse = "-")
}

# # Now summarise blocks within species_comparisons
# ggplot(raw_collinearity_groups[raw_collinearity_groups$group_length < 100,],aes(species_pair,y=PASS_genepairs))+
#   geom_boxplot()

# Mark groups that pass/fail quality checks
raw_collinearity_groups$pass_fail <- "FAIL"
raw_collinearity_groups$pass_notpass_ratio <- raw_collinearity_groups$PASS_genepairs/raw_collinearity_groups$NOTPASS_genepairs
raw_collinearity_groups[raw_collinearity_groups$PASS_genepairs >= 2 & raw_collinearity_groups$PASS_genepairs > raw_collinearity_groups$NOTPASS_genepairs,"pass_fail"] <- "PASS"
for(species in unique(raw_collinearity_groups$species_pair)){
  print(species)
  print(nrow(raw_collinearity_groups[raw_collinearity_groups$species_pair == species & raw_collinearity_groups$pass_fail == "PASS",])/nrow(raw_collinearity_groups[raw_collinearity_groups$species_pair == species,]))
}

# Now summarise blocks within species_comparisons
ggplot(raw_collinearity_groups,aes(species_pair,y=log(pass_notpass_ratio)))+
  geom_boxplot()

# Assess pSONIC tethers ---------------------------------------------------
tethers <- read.table(paste0(results_dir,"pSONIC.TetherSetsFromOrthoFinder.csv"),sep="\t")[,1:4]
colnames(tethers) <- colnames(OF_groups)[2:ncol(OF_groups)]

# Count the number of representa
tether_counts <- count_number_genes_per_set(tethers)
to_plot <- data.frame(tether_count = colSums(tether_counts),
                      species = colnames(tethers))
ggplot(to_plot,aes(x=species,y=tether_count))+
  geom_bar(stat="identity")+
  ggtitle("Genes included in MCScanX Tethers")

# Visualise pSONIC summaries ----------------------------------------------


