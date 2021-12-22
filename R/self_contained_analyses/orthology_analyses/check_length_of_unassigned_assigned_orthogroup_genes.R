# Analysis of whether filtering orthofinder genes prior to OF2 is pointless
lib <- c("data.table","ggplot2","viridis","ggridges","dplyr","pbmcapply","ggridges","cowplot")
sapply(lib,library,character.only=T)

# Orthofinder2 res...
OF_res <- "outputs/orthology/Results_210706_orthofinder_test/"

# Get all the proteomes...
prot_fai <- lapply(paste0(OF_res,"/",list.files(OF_res,pattern = "fai")),read.table,header=F)
names(prot_fai) <- gsub(".faa.fai","",list.files(OF_res,pattern = "fai"))

# Fetch all the orthogroups
OG <- read.table(paste0(OF_res,"/Orthogroups/Orthogroups.tsv"),sep="\t",header=T)
unassigned <- read.table(paste0(OF_res,"/Orthogroups/Orthogroups_UnassignedGenes.tsv"),sep="\t",header=T)
  
# Get all the sequence IDs
seqID <- read.table(paste0(OF_res,"/SequenceIDs.txt"),fill=T)

# Work out our species codes
all_species <- list.files(OF_res,pattern = "fai")
all_species <- gsub(".faa.fai","",all_species)
all_species <- all_species[all_species %in% colnames(OG)]
prot_fai <- prot_fai[all_species]

# For each species, we now want to collect the sizes of unassigned genes...
species_unassigned_lengths <- lapply(all_species,function(species){
  
  # Fetch the unassigned
  message(paste0(">>> STARTING ",species))
  tmp_unassigned <- na.omit(unassigned[,species])
  tmp_unassigned <- gsub("transcript_","transcript:",tmp_unassigned)
  tmp_unassigned <- tmp_unassigned[tmp_unassigned != ""]
  
  # Find the length of these in the relevant fai
  all_matches <- pbmclapply(tmp_unassigned,function(gene){
   # tmp_length <- prot_fai[[species]][grep(gene,prot_fai[[species]]$V1),"V2"]
   # if(length(tmp_length) == 1){
   #   return(tmp_length)
   # } else {
   #   return(NA)
   # }
    tmp_match <- grep(gene,prot_fai[[species]]$V1)
    if(length(tmp_match) == 1){
      return(tmp_match)
    } else {
      return(NA)
    }
  },mc.cores = 6)
  
  # Can now just fetch all lengths and also remove the unassigned from full prot
  assigned_prot <- prot_fai[[species]][-na.omit(unlist(all_matches)),]
  all_lengths <- prot_fai[[species]][unlist(all_matches),"V2"]
  
  # Plot alongside the original proteome...
  to_plot <- data.frame(length=c(all_lengths,assigned_prot$V2),
                        group=c(rep("Unassigned",length(all_lengths)),rep("Assigned",length(assigned_prot$V2))),
                        species=species)
  g1 <- ggplot(to_plot,aes(x=length,y=group))+
   # geom_density(alpha=0.5)+
    stat_density_ridges(quantile_lines = TRUE, quantiles = 2)+
    xlim(0,500)+
    ggtitle(species)+
    labs(y="",x="Protein length (AA)")
  
  return(list(g1,to_plot))
})

# All figs
unassigned_figs <- lapply(species_unassigned_lengths,'[[',1)
plot_grid(plotlist=unassigned_figs)

# Fetch all the plotting dfs
all_species_lengths <- lapply(species_unassigned_lengths,'[[',2)
for(x in all_species_lengths){
  # Get median
  tmp_median <- median(x[x$group == "Unassigned","length"],na.rm = T)
  
  # What prop of full is lower than this...
  lower_than_median <- nrow(x[x$group == "Assigned" & x$length < tmp_median,])/nrow(x[x$group == "Assigned",])
  lower_than_50 <- nrow(x[x$group == "Assigned" & x$length < 50,])/nrow(x[x$group == "Assigned",])
  lower_than_50_un <- nrow(x[x$group == "Unassigned" & x$length < 50,])/nrow(x[x$group == "Unassigned",])
  
  # Print res
  print(unique(x$species))
  print(paste0(round(lower_than_median,3)*100,"% Assigned genes lower than Unassigned median of ",tmp_median))
  print(paste0(round(lower_than_50,3)*100,"% Assigned genes lower than AA filter of 50"))
  print(paste0(round(lower_than_50_un,3)*100,"% Unassigned genes lower than AA filter of 50"))
  
}

