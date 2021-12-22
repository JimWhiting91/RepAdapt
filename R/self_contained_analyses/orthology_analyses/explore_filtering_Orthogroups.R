# Look at lop-sided orthogroups
lib <- c("ggplot2","data.table","moments","tidyverse")
lapply(lib,library,character.only=T)

# Fetch gene counts
genecounts <- read.table("outputs/Orthogroups.GeneCount.tsv",header=T)
genecounts <- genecounts[,colnames(genecounts) != "Total"]
rownames(genecounts) <- genecounts$Orthogroup
genecounts_clean <- genecounts[,colnames(genecounts) != "Orthogroup"]

# Calculate per-row skew
genecounts$ortho_skew <- apply(genecounts_clean,1,skewness)

# Calculate max distance from the median
max_median_distance <- function(vec){
  med <- median(vec)
  max_dist <- max(abs(vec-med))
  max_dist
}
genecounts$max_median_distance <- apply(genecounts_clean,1,max_median_distance)

# Associate these
ggplot(genecounts,aes(max_median_distance,ortho_skew))+
  geom_point()+
  geom_vline(xintercept = 10,colour="red2")

# Associate these after cutting
ggplot(genecounts[genecounts$max_median_distance < 10,],aes(max_median_distance,ortho_skew))+
  geom_point()


# Filtering thresholds and OG retention -----------------------------------
filter_mat <- matrix(ncol = 2,nrow=20)
filter_mat[,1] <- 20:1

for(i in 1:nrow(filter_mat)){
  # How many orthogroups are left if we filter for any case in which an orthogroup includes more than N genes per species...
  genecounts_filtered <- genecounts_clean
  genecounts_filtered[genecounts_filtered > filter_mat[i,1]] <- NA
  genecounts_filtered <- na.omit(genecounts_filtered)
  
  # What percentage do we retain
  filter_mat[i,2] <- round(nrow(genecounts_filtered)/nrow(genecounts_clean)*100,3)
}

as.data.frame(filter_mat) %>%
  ggplot(aes(V1,V2))+
  geom_line()+
  labs(y="% OG retained",x="Filtering threshold (remove OG where any species has > x)")
