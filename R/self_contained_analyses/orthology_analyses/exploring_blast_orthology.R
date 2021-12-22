# Pairwise orthology based on OF2 Blast results...
lib <- c("lostruct","data.table","ggplot2","viridis","ggridges","dplyr","patchwork","readr")
sapply(lib,library,character.only=T)

# Function Library --------------------------------------------------------

# Where are the blast results?
blast_dir <- "outputs/orthology/Results_210913_17_genomes_with_conifers_noAA_filter/WorkingDirectory/"
blast_files <- list.files(blast_dir,pattern="Blast")

# Read one in to explore in-depth
blast_res <- data.frame(read_table(paste0(blast_dir,"/",blast_files[2]),progress = T,col_names = F))
colnames(blast_res) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","start","ssend","evalue","bitscore")

# Let's also just make up some random WZA scores for these...
wza1 <- data.frame(wza=rnorm(length(unique(blast_res$qseqid)),mean=0,sd=3),
                   qseqid=unique(blast_res$qseqid))
wza2 <- data.frame(wza=rnorm(length(unique(blast_res$sseqid)),mean=0,sd=3),
                   sseqid=unique(blast_res$sseqid))

# And get their ranks...
wza1$rank <- rank(wza1$wza)
wza2$rank <- rank(wza2$wza)

# Transform the evalues to allow for min values...
blast_res$evalue_trans <- -log10(blast_res$evalue + min(blast_res$evalue[blast_res$evalue > 0]))

# Let's explore this
ggplot(blast_res[1:1000,],aes(bitscore,evalue_trans))+
  geom_point(alpha=0.5)+
  geom_density2d()

# Count occurrence of query seqs...
hist(table(blast_res$qseqid))
hist(table(blast_res$sseqid))
sum(table(blast_res$qseqid))
sum(table(blast_res$sseqid))

# How many unique of each...
length(unique(blast_res$qseqid))
length(unique(blast_res$sseqid))

# View the distribution of blast hits for the most common sseq
most_common_sseq <- names(which(table(blast_res$sseqid) == max(table(blast_res$sseqid))))
hist(blast_res[blast_res$sseqid == most_common_sseq,"bitscore"])
ggplot(blast_res[blast_res$sseqid == most_common_sseq,],aes(bitscore,evalue_trans))+
  geom_point(alpha=0.5)+
  geom_density2d()

# Playing about with matrices
ggplot(blast_res[1:1000,],aes(qseqid,sseqid,fill=evalue_trans))+
  geom_tile()+
  scale_fill_viridis(option="A")

# Let's try and quantify some aspects of different genes...
qseq_sums <- pbmcapply::pbmclapply(unique(blast_res$qseqid[1:1000]),function(id){
  
  # Subset
  tmp <- blast_res[blast_res$qseqid == id, ]
  
  # For each of 
  
},mc.cores=4)

# # Convert to wide matrix and quantify sparcity...
# test <- data.frame(tidyr::pivot_wider(blast_res[,c("qseqid","sseqid","evalue_trans")],names_from="qseqid",values_from="evalue_trans"))
# sparcity <- sum(is.na(test))/(ncol(test)*nrow(test))


# Attempt at "correlation of repeatability" approach ----------------------
blast_res_clean <- blast_res[,c("qseqid","sseqid","bitscore")]

# And we want the ranking of every gene in both seq columns...
blast_res_clean$qseq_ran
blast_res_clean <- merge(blast_res_clean,wza1[,c("qseqid","wza")],"qseqid")
colnames(blast_res_clean)[ncol(blast_res_clean)] <- "q_wza"
blast_res_clean <- merge(blast_res_clean,wza2[,c("sseqid","wza")],"sseqid")
colnames(blast_res_clean)[ncol(blast_res_clean)] <- "s_wza"

# Now perform a weighted correlation...
library(wCorr)
test_corr <- weightedCorr(blast_res_clean$q_wza,blast_res_clean$s_wza,method = "Pearson",weights = blast_res_clean$bitscore)


# Attempt at Z-score based approach ---------------------------------------
library(pbmcapply)
Z_perms = 1000

wza1_outliers <- wza1[wza1$wza > quantile(wza1$wza,0.99),]
wza2_outliers <- wza2[wza2$wza > quantile(wza2$wza,0.99),]

# It might be quicker to just expand out wza1_outliers once rather than repeatedly sample within the loops...
blast_res_expanded <- data.frame(rbindlist(lapply(unique(wza1_outliers$qseqid),function(qseqid){
  tmp <- blast_res_clean[blast_res_clean$qseqid==qseqid,]
  tmp[rep(seq_len(nrow(tmp)), times = tmp$bitscore), ]
})))

perm_res <- pbmclapply(1:Z_perms,function(iter){
  set.seed(iter)
  # # Based on bitscores, randomly draw 1 gene from species 2 based on the outliers from species 1
  # tmp_ortho_fetch <- sapply(wza1_outliers$qseqid,function(qseqid){
  #   sample(blast_res_clean[blast_res_clean$qseqid==qseqid,"sseqid"],prob = blast_res_clean[blast_res_clean$qseqid==qseqid,"bitscore"],1)
  # })
  
  # Take a sample from the expanded blast_res_clean
  sample_tmp <- blast_res_expanded[sample(1:nrow(blast_res_expanded)),]
  sample_tmp <- sample_tmp[!(duplicated(sample_tmp$qseqid)),]
  tmp_ortho_fetch <- sample_tmp$sseqid
  
  # How many of these are similarly in the top 1% - We expect 1%...
  length(tmp_ortho_fetch[tmp_ortho_fetch %in% wza2_outliers$sseqid])/nrow(wza2_outliers)
},mc.cores=6)

# Our Z score is then something along the lines of the median of this distribution compared with the null or something else?
hist(unlist(perm_res))
median(unlist(perm_res))/0.01

