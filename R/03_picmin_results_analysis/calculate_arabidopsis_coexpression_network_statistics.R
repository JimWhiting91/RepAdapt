# Calculate node-level statistics for Arabdidopsis coexpression network...
lib = c("qvalue","ggplot2","pbmcapply","biomaRt","org.At.tair.db","ExpressionAtlas","S4Vectors", "IRanges", "GenomicRanges", "SummarizedExperiment","readr","data.table")
sapply(lib,library,character.only=T)
n_cores = 6

# What run?
run_name = "220927"
output_name = "COMBINED25species_updatedOF_221213"

# Fetch the OG map Athal
OG_map_Athal = readRDS(paste0("data/OG_map_Athal_",run_name,".rds"))

# First we make a table of TAIR gene ID vs Entrez through biomart
# Set up Ensembl
ensembl <- useMart(biomart = "plants_mart",host="https://plants.ensembl.org")
athal_biomart <- useDataset("athaliana_eg_gene",mart=ensembl)

# Fetch whole biomart entry for thaliana
athal_universe <- getBM(attributes = c("tair_locus","entrezgene_id"),
                        mart=athal_biomart)
colnames(athal_universe) <- c("TAIR_gene","ENTREZ_gene")

# Merge these into our map
OG_map_Athal_entrez = data.table(merge(OG_map_Athal,athal_universe,by = "TAIR_gene"))

# Read in all of the coexpression files and bind together
coexpression_Zcutoff = 2.333
neg_coexpression_Zcutoff = -5
all_coexpression_files = list.files("data/coexpression_data/")
all_coexpression_files = grep("zip",all_coexpression_files,invert = T,value = T)
gene_coexpression = pbmclapply(all_coexpression_files,function(input){
  
  # print(input)
  
  # Fetch the file
  tmp = fread(paste0("data/coexpression_data/",input))
  colnames(tmp) = c("ENTREZ_gene","CoEx_Z")
  
  # Merge and reduce down to 'signif' coexpression
  tmp2 = merge(tmp[CoEx_Z >= coexpression_Zcutoff | 
                     CoEx_Z <= neg_coexpression_Zcutoff,],OG_map_Athal_entrez[,c("TAIR_gene","ENTREZ_gene")])
  focal_TAIR = OG_map_Athal_entrez[ENTREZ_gene == input,TAIR_gene]
  if(length(focal_TAIR) == 1){
    tmp2$focal_TAIR = focal_TAIR
    return(tmp2)
  } else {
    return(NULL)
  }
},mc.cores = n_cores)

# In the darkness, bind them
all_gene_coexpression = rbindlist(gene_coexpression)

# Convert to matrix and then network...
coex_matrix = as.matrix(all_gene_coexpression[,.(focal_TAIR,TAIR_gene)])
coex_ig <- graph.edgelist(coex_matrix , directed=TRUE)
coex_ig2 <- graph.data.frame(data.table(from = all_gene_coexpression$focal_TAIR,
                                        to = all_gene_coexpression$TAIR_gene,
                                        weight = abs(all_gene_coexpression$CoEx_Z)))

# Fetch stats
# Betweenness, cutoff of -1 ensures that no cutoff is used for betweeness
coex_betweenness = estimate_betweenness(coex_ig2,cutoff = -1)
coex_betweenness = data.table(TAIR_gene = names(coex_betweenness),
                              node_betweenness = coex_betweenness)
# Node degree and strength
coex_strength = all_gene_coexpression[,.(node_strength = sum(abs(CoEx_Z)),
                                         node_degree = nrow(.SD)),by = focal_TAIR]
colnames(coex_strength)[1] = c("TAIR_gene")
# Closeness
coex_closeness = closeness(coex_ig2)
coex_closeness = data.table(TAIR_gene = names(coex_closeness),
                              node_closeness = coex_closeness)

# Merge all of these...
coex_merge = merge(merge(coex_betweenness,coex_strength),coex_closeness)

# Save this
saveRDS(coex_merge,
        paste0("outputs/",run_name,"_Athal_coexpression_node_stats.rds"))


# Repeat for Medicago... --------------------------------------------------
OG_map_Mtrunc = data.table(readRDS(paste0("data/OG_map_Mtrunc_",output_name,".rds")))


# First we make a table of TAIR gene ID vs Entrez through biomart
# Set up Ensembl
ensembl <- useMart(biomart = "plants_mart",host="https://plants.ensembl.org")
Mtrunc_biomart <- useDataset("mtruncatula_eg_gene",mart=ensembl)

# Fetch whole biomart entry for thaliana
Mtrunc_universe <- getBM(attributes = c("ensembl_gene_id","entrezgene_id"),
                        mart=Mtrunc_biomart)
colnames(Mtrunc_universe) <- c("biomart_gene","ENTREZ_gene")

# Merge these into our map
OG_map_Mtrunc_entrez = data.table(merge(OG_map_Mtrunc[,c("biomart_gene","Orthogroup")],Mtrunc_universe,by = "biomart_gene"))

# Read in all of the coexpression files and bind together
coexpression_Zcutoff = 2.333
neg_coexpression_Zcutoff = -5
all_coexpression_files = list.files("data/medicago_coexpression_data/")
all_coexpression_files = grep("zip",all_coexpression_files,invert = T,value = T)
gene_coexpression = pbmclapply(all_coexpression_files,function(input){
  
  # print(input)
  
  # Fetch the file
  tmp = fread(paste0("data/medicago_coexpression_data/",input))
  colnames(tmp) = c("ENTREZ_gene","CoEx_Z")
  
  # Merge and reduce down to 'signif' coexpression
  tmp2 = merge(tmp[CoEx_Z >= coexpression_Zcutoff | 
                     CoEx_Z <= neg_coexpression_Zcutoff,],OG_map_Mtrunc_entrez[,c("biomart_gene","ENTREZ_gene")])
  focal_biomart = OG_map_Mtrunc_entrez[ENTREZ_gene == input,biomart_gene]
  if(length(focal_biomart) == 1){
    tmp2$focal_biomart = focal_biomart
    return(tmp2)
  } else {
    return(NULL)
  }
},mc.cores = n_cores)

# In the darkness, bind them
all_gene_coexpression = rbindlist(gene_coexpression)

# Convert to matrix and then network...
library(igraph)
coex_matrix = as.matrix(all_gene_coexpression[,.(focal_biomart,biomart_gene)])
coex_ig <- graph.edgelist(coex_matrix , directed=TRUE)
coex_ig2 <- graph.data.frame(data.table(from = all_gene_coexpression$focal_biomart,
                                        to = all_gene_coexpression$biomart_gene,
                                        weight = abs(all_gene_coexpression$CoEx_Z)))

# Fetch stats
# Betweenness, cutoff of -1 ensures that no cutoff is used for betweeness
coex_betweenness = estimate_betweenness(coex_ig2,cutoff = -1)
coex_betweenness = data.table(biomart_gene = names(coex_betweenness),
                              node_betweenness = coex_betweenness)
# Node degree and strength
coex_strength = all_gene_coexpression[,.(node_strength = sum(abs(CoEx_Z)),
                                         node_degree = nrow(.SD)),by = focal_biomart]
colnames(coex_strength)[1] = c("biomart_gene")
# Closeness
coex_closeness = closeness(coex_ig2)
coex_closeness = data.table(biomart_gene = names(coex_closeness),
                            node_closeness = coex_closeness)

# Merge all of these...
coex_merge = merge(merge(coex_betweenness,coex_strength),coex_closeness)

# Save this
saveRDS(coex_merge,
        paste0("outputs/",run_name,"_Mtrunc_coexpression_node_stats.rds"))

