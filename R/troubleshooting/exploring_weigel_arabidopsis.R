# Exploratory analysis of weigel Athaliana SNPs..
library(vcfR)
library(data.table)
library(ggplot2)
library(dplyr)
library(cowplot)
library(regioneR)

# Fetch subset
vcf <- read.vcfR("data/VCFs/Athaliana_weigel_Athal_NC_003070.9_all_inds.vcf.gz")

# Recreate fucked up distribution - yep
vcf_pos <- data.frame(chr="Athal_NC_003070.9",
                      pos=as.integer(vcf@fix[,2]))
hist(vcf_pos$pos)

# Plot distributions of snp scores...
GQ_mat <- extract.gt(vcf,element = "GQ")
hist(GQ_mat)

AD_mat <- extract.gt(vcf,element = "AD")
AD_mat_sub <- AD_mat[1:50000,]
for(i in 3:ncol(AD_mat)){
  print(paste0("Starting col:",i))
  firstAD <- as.integer(sapply(strsplit(na.omit(AD_mat_sub[,i]),","),'[[',1))
  secondAD <- as.integer(sapply(strsplit(na.omit(AD_mat_sub[,i]),","),'[[',2))
  AD_mat_sub[!(is.na(AD_mat_sub[,i])),i] <- firstAD+secondAD
}
hist(as.numeric(AD_mat_sub)[as.numeric(AD_mat_sub) < 20])
median(as.numeric(AD_mat_sub),na.rm = T)
mean(as.numeric(AD_mat_sub),na.rm = T)

# What proportion of SNPs are <5 AD
table(as.numeric(AD_mat_sub) < 5)[2]/sum(table(as.numeric(AD_mat_sub) < 5))

# Find prop of missing data
table(is.na(GQ_mat))[2]/sum(table(is.na(GQ_mat)))

# There's lots of missing data, what's the per-site rate of missing
GT_mat <- extract.gt(vcf,element = "GT")
persite_missing <- sapply(1:nrow(GT_mat),function(row){
  sum(is.na(GT_mat[row,]))/ncol(GT_mat)
})

# Plot this
plot_missing <- data.frame(pos=as.integer(vcf@fix[,2]),
                           missing=persite_missing)
ggplot(plot_missing,aes(pos,missing))+
  geom_point()+
  geom_smooth()

# Average into windows...
winds <- seq(0,max(plot_missing$pos),1000)
plot_missing_winds <- data.frame(rbindlist(lapply(winds,function(wind){
  data.frame(start=wind,
             missing_avg=median(plot_missing[plot_missing$pos >= wind &
                                               plot_missing$pos < wind+1000,"missing"]))
})))

plot_missing_winds %>%
  ggplot(aes(start,missing_avg))+
  geom_point()+
  stat_density_2d()



# Compare with Ahalleri SNPs ----------------------------------------------
ahalleri_snps <- data.frame(fread("data/VCFs/rellstab_Ahal_pooled-varscan_all_bedfiles_SNP_filtered_NoMiss0.25_maf0.05.txt"))
ahalleri_snps <- ahalleri_snps[,1:2]
ahalleri_snps$rand <- rnorm(nrow(ahalleri_snps),0,1)

# Plot histogram over biggest chr
chr_max <- data.frame(ahalleri_snps %>% group_by(CHROM) %>% summarise(max_pos=max(POS)))
chr_max <- chr_max[order(-chr_max$max_pos),]

hist(ahalleri_snps[ahalleri_snps$CHROM == "Ahall_FJVB01000011.1","POS"])
ahalleri_snps[ahalleri_snps$CHROM == "Ahall_FJVB01000011.1",] %>%
  ggplot(aes(POS,rand))+geom_point()

# Make bedfile of all genes where there is no
Ahalleri_prot <- data.frame(fread("data/reference_genomes/A.halleri_gemmifera_reference/Ahalleri_proteome_to_OF_id_map.txt"))
Ahalleri_prot$gea_gene <- paste0(Ahalleri_prot$seqid,":",Ahalleri_prot$start,"-",Ahalleri_prot$end)
Ahalleri_prot_regions <- toGRanges(Ahalleri_prot[,c("seqid","start","end")])
ahalleri_snps$POS2 <- ahalleri_snps$POS+1
ahalleri_snps_regions <- toGRanges(ahalleri_snps[,c("CHROM","POS","POS2")])
overlap <- overlapRegions(Ahalleri_prot_regions,ahalleri_snps_regions)
overlap$gene <- paste0(overlap$chr,":",overlap$startA,"-",overlap$endA)
length(unique(overlap$gene))/nrow(Ahalleri_prot)

# Which genes are missing from SNPs...
missing_genes <- Ahalleri_prot$gea_gene[!(Ahalleri_prot$gea_gene %in% unique(overlap$gene))]
missing_genes_bed <-  tidyr::separate(data.frame(missing_genes),col = "missing_genes",into=c("chr","start","end"),sep = ":|-")
write.table(missing_genes_bed,
            "outputs/rellstab_Ahalleri_missing_genes.bed",
            row.names = F,col.names = F,sep="\t",quote = F)

# Compare with original Weigel 1,135 genome paper SNPs --------------------
weigel_snps <- data.frame(fread("data/VCFs/1001_genomes_chr1_snps.pos"))
weigel_snps$rand <- rnorm(nrow(weigel_snps),0,1)

# Compare hists of snp density
library(cowplot)
plot_grid(
  ggplot(weigel_snps,aes(x=V2))+geom_histogram(bins=100),
  ggplot(vcf_pos,aes(x=pos))+geom_histogram(bins=100),
  ncol=1)

# Compare weigel snps to Athaliana genome genes
Athaliana_prot <- data.frame(fread("data/reference_genomes/A.thaliana_reference/Athaliana_proteome_to_OF_id_map.txt"))
Athaliana_prot_regions <- toGRanges(Athaliana_prot[Athaliana_prot$seqname=="Athal_NC_003070.9",c("seqname","start","end")])
weigel_snps$pos2 <- weigel_snps$V2+1
weigel_snps$chr <- "Athal_NC_003070.9"
weigel_snps_regions <- toGRanges(weigel_snps[,c("chr","V2","pos2")])
overlap <- overlapRegions(Athaliana_prot_regions,weigel_snps_regions)
overlap$gene <- paste0(overlap$chr,":",overlap$startA,"-",overlap$endA)
length(unique(overlap$gene))/nrow(Athaliana_prot[Athaliana_prot$seqname=="Athal_NC_003070.9",])


# Fetch all the accession data --------------------------------------------
weigel_accessions <- read.csv("metadata/dataset_SRA_ENA/weigel_Athaliana_all_accession_info",header=F)

# Count up countries
colnames(weigel_accessions) <- c("Accession_ID","Sequenced","Name","Country","Location?","Lat","Long","Collector","","CS_Number","Admixture_Group","","plot_code")
country_counts <- table(weigel_accessions$Country)
sort(country_counts)

# Filter the swedish data from teh sra data and export the download codes...
# Get swedish
swedish_samples <- weigel_accessions[weigel_accessions$Country=="SWE","Accession_ID"]

# Get Iberian, but exclude "relicts"
iberian_samples <- weigel_accessions[weigel_accessions$Country %in% c("POR","ESP") &
                                       !(weigel_accessions$Admixture_Group %in% c("relict")),"Accession_ID"]

# Get Iberian, but exclude "relicts"
usa_samples <- weigel_accessions[weigel_accessions$Country %in% c("USA") &
                                       !(weigel_accessions$Admixture_Group %in% c("relict")),"Accession_ID"]

# Combine country and admix group...
weigel_accessions$country_admix <- paste0(weigel_accessions$Country,"-",weigel_accessions$Admixture_Group)
country_admix_counts <- table(weigel_accessions$country_admix)
sort(country_admix_counts)

sra_data <- read.csv("metadata/dataset_SRA_ENA/weigel_Athalian_SRR_run_info.csv")
swedish_sra <- sra_data[sra_data$SampleName %in% swedish_samples,"Run"]
write.table(data.frame(swedish_sra),
            "metadata/dataset_SRA_ENA/weigel_Athaliana_swedish_accession_sra_download_list.txt",
            row.names = F,col.names = F,quote = F)

iberia_sra <- sra_data[sra_data$SampleName %in% iberian_samples,"Run"]
write.table(data.frame(iberia_sra),
            "metadata/dataset_SRA_ENA/weigel_Athaliana_iberia_accession_sra_download_list.txt",
            row.names = F,col.names = F,quote = F)

usa_sra <- sra_data[sra_data$SampleName %in% usa_samples,"Run"]
write.table(data.frame(usa_sra),
            "metadata/dataset_SRA_ENA/weigel_Athaliana_USA_accession_sra_download_list.txt",
            row.names = F,col.names = F,quote = F)

# Count up admixture groups...
admix_counts <- table(weigel_accessions$Admixture_Group)
sort(admix_counts)


# Going through BAM file info ---------------------------------------------
# Swedish Arabidopsis
bamstats <- read.table("metadata/snpcalling/weigel_Athaliana_final_bam_stats.txt",comment.char = "#",header = T,row.names = NULL)
colnames(bamstats) <- c(colnames(bamstats)[2:ncol(bamstats)],"individual")
bamstats$dataset <- "weigel_Athaliana_SWE"

# Plot out coverages...
bamstats %>%
  ggplot(aes(x=MEAN_COVERAGE,fill=CATEGORY))+
  geom_histogram(bins=100)+
  facet_wrap(~CATEGORY,ncol=1)
bamstats %>%
  ggplot(aes(x=MEDIAN_COVERAGE,fill=CATEGORY))+
  geom_histogram(bins=100)+
  facet_wrap(~CATEGORY,ncol=1)

# What is the size of non-zero regions...
nonzero_prop <- sapply(unique(bamstats$individual),function(ind) {
  bamstats[bamstats$individual==ind & bamstats$CATEGORY=="NON_ZERO_REGIONS","GENOME_TERRITORY"] /
    bamstats[bamstats$individual==ind & bamstats$CATEGORY=="WHOLE_GENOME","GENOME_TERRITORY"]  
})
hist(nonzero_prop)

# Various PCT metrics
PCT_metrics <- grep("PCT_EX",colnames(bamstats),value=T)
PCT_figs <- lapply(PCT_metrics,function(PCT){
bamstats[bamstats$CATEGORY=="NON_ZERO_REGIONS",] %>%
  ggplot(aes(x=bamstats[bamstats$CATEGORY=="NON_ZERO_REGIONS",PCT]))+
  geom_histogram(bins=100)+
  labs(x=PCT)
})
cowplot::plot_grid(plotlist=PCT_figs,ncol=2)
  
# Mtrunc
bamstats2 <- read.table("metadata/snpcalling/rellstab_pool_v1_bamstats.txt",comment.char = "#",header = T,row.names = NULL)
colnames(bamstats2) <- c(colnames(bamstats2)[2:ncol(bamstats2)],"individual")
bamstats2$dataset <- "tiffin_Mtrunc"

# Plot out coverages...
bamstats2 %>%
  ggplot(aes(x=MEAN_COVERAGE,fill=CATEGORY))+
  geom_histogram(bins=100)+
  facet_wrap(~CATEGORY,ncol=1)
bamstats2 %>%
  ggplot(aes(x=MEDIAN_COVERAGE,fill=CATEGORY))+
  geom_histogram(bins=100)+
  facet_wrap(~CATEGORY,ncol=1)

# What is the size of non-zero regions...
nonzero_prop <- sapply(unique(bamstats2$individual),function(ind) {
  bamstats2[bamstats2$individual==ind & bamstats2$CATEGORY=="NON_ZERO_REGIONS","GENOME_TERRITORY"] /
    bamstats2[bamstats2$individual==ind & bamstats2$CATEGORY=="WHOLE_GENOME","GENOME_TERRITORY"]  
})
hist(nonzero_prop)

# Various PCT metrics
PCT_metrics <- grep("PCT_EX",colnames(bamstats2),value=T)
PCT_figs2 <- lapply(PCT_metrics,function(PCT){
  bamstats2[bamstats2$CATEGORY=="NON_ZERO_REGIONS",] %>%
    ggplot(aes(x=bamstats2[bamstats2$CATEGORY=="NON_ZERO_REGIONS",PCT]))+
    geom_histogram(bins=100)+
    labs(x=PCT)
})
cowplot::plot_grid(plotlist=PCT_figs2,ncol=2)


# Recalled SNPs sweden only... --------------------------------------------
# Compare weigel snps to Athaliana genome genes
weigel_swe_snps <- data.frame(fread("data/VCFs/weigel_Athaliana_SWE.snps"))
Athaliana_prot <- data.frame(fread("data/reference_genomes/A.thaliana_reference/Athaliana_proteome_to_OF_id_map.txt"))
Athaliana_prot_regions <- toGRanges(Athaliana_prot[,c("seqname","start","end")])
weigel_swe_snps$pos2 <- weigel_swe_snps$V2+1
weigel_swe_snps_regions <- toGRanges(weigel_swe_snps[,c("V1","V2","pos2")])
overlap <- overlapRegions(Athaliana_prot_regions,weigel_swe_snps_regions)
overlap$gene <- paste0(overlap$chr,":",overlap$startA,"-",overlap$endA)
length(unique(overlap$gene))/nrow(Athaliana_prot)

hist(weigel_swe_snps[weigel_swe_snps$V1 == "Athal_NC_003070.9","V2"])


