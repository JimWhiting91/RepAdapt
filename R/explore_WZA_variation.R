#######################################################################
# Quick script to explore variation in WZA results across variables

# Set up environment
lib <- c("parallel","cowplot","data.table","ggplot2","ggExtra")
lapply(lib,library,character.only=T)

# Fetch WZA results
res <- list.files("outputs/GEA_res/Arabidopsis_halleri_Kubota_Individual/",pattern = "WZA")

# Set up climate vars
climate_vars <- c("mean_temp","mean_diurnal","isothermality","temp_seasonality","max_temp_warmest_month",
                  "min_temp_coldest_month","temp_range","mean_temp_wet_quarter","mean_temp_dry_quarter","mean_temp_warm_quarter",
                  "mean_temp_cold_quarter","annual_precip","precip_wet_month","precip_dry_month","precip_seasonality","precip_wet_quarter",
                  "precip_dry_quarter","precip_warm_quarter","precip_cold_quarter")

# Read in all the res
wza_res <- data.frame(rbindlist(lapply(res,function(x){
  fread(paste0("outputs/GEA_res/Arabidopsis_halleri_Kubota_Individual/",x))
})))

# Visualise distributions
ggplot(wza_res,aes(x=log10(weiZ),y=climate_var))+
  geom_density_ridges()

# Cor mat
res_mat <- matrix(ncol=length(climate_vars),nrow=nrow(wza_res)/length(climate_vars))
for(i in 1:ncol(res_mat)){
  res_mat[,i] <- wza_res[wza_res$climate_var == climate_vars[i],"weiZ"]
}
res_corr <- cor(res_mat)
rownames(res_corr) <- climate_vars
colnames(res_corr) <- climate_vars

library(corrplot)
corrplot(res_corr,method = "ellipse",order="hclust")

# For each climate var. what's the maximum weighted Z
for(climate_var in climate_vars){
  print(max(wza_res[wza_res$climate_var==climate_var,"weiZ"]))
}

# For each climate var. what's the variance of the distribution
for(climate_var in climate_vars){
  print((wza_res[wza_res$climate_var==climate_var,"weiZ"]))
}

#########################################################################################################
##### Compare WZA and TC approaches ######
wza_res <- list.files("outputs/GEA_res/Arabidopsis_halleri_Kubota_Individual/",pattern = "WZA")
TC_res <- list.files("outputs/GEA_res/Arabidopsis_halleri_Kubota_Individual/",pattern = "TC")

# Fetch the gene info with pbar_qbar
gene_summary <- data.frame(fread("outputs/GEA_res/Arabidopsis_halleri_Kubota_Individual/GEA_snpstats_pergene.tsv"))

# Do for all climate vars
climate_TC_WZA <- lapply(climate_vars,function(var){
  print(var)
  
  # Compare the annual both_geaitation results
  both_gea <- cbind(data.frame(fread(paste0("outputs/GEA_res/Arabidopsis_halleri_Kubota_Individual/",var,"_WZA_pergene.tsv"))),
                    data.frame(fread(paste0("outputs/GEA_res/Arabidopsis_halleri_Kubota_Individual/",var,"_TC_pergene.tsv"))))
  
  # Mark that those agree and those don't
  wza_cutoff <- 4
  both_gea$overlapping <- NA
  both_gea[both_gea$TC_score > 0,"overlapping"] <- "TC"
  both_gea[both_gea$weiZ > wza_cutoff,"overlapping"] <- "WZA"
  both_gea[both_gea$TC_score > 0 & both_gea$weiZ > wza_cutoff,"overlapping"] <- "Both"
  # Get overlap
  overlap_stats <- table(both_gea$overlapping)
  overlap_title <- paste0(var," overlap = WZA:(",overlap_stats["Both"],"/",sum(overlap_stats[c("Both","WZA")]),")",
                          " TC:(",overlap_stats["Both"],"/",sum(overlap_stats[c("Both","TC")]),")")
  
  both_gea <- both_gea[,-c(1,3)]
  weiZ_TC <- ggplot(both_gea,aes(weiZ,TC_score,colour=overlapping))+
    geom_point()+
    geom_hline(yintercept=0,colour="red2")+
    geom_vline(xintercept=wza_cutoff,colour="red2")+
    ggtitle(overlap_title)
  
  # Return ovelrap
  overlap_out <- data.frame(overlap_stats,
                            var)
  
  # Plot these relative to corr
  weiz_corr <- ggplot(both_gea[!(is.na(both_gea$overlapping)),],aes(weiZ,mean_corr,colour=overlapping))+
    geom_point()+
    geom_vline(xintercept=wza_cutoff,colour="red2")
  weiz_corr_marginal <- ggMarginal(weiz_corr, groupColour = TRUE, groupFill = TRUE)
  
  TC_corr <- ggplot(both_gea[!(is.na(both_gea$overlapping)),],aes(TC_score,mean_corr,colour=overlapping))+
    geom_point()+
    geom_vline(xintercept=0,colour="red2")
  TC_corr_marginal <- ggMarginal(TC_corr, groupColour = TRUE, groupFill = TRUE)
  
  # Merge with pbar_qbar snp_count
  both_gea_merge <- merge(both_gea,gene_summary,by="gene_id")
  
  # Plot these relative to pbar_qbar
  weiz_pbar_qbar <- ggplot(both_gea_merge[!(is.na(both_gea_merge$overlapping)),],aes(mean_pbar_qbar,weiZ,colour=overlapping))+
    geom_point()+
    geom_hline(yintercept=wza_cutoff,colour="red2")
  weiz_pbar_qbar_marginal <- ggMarginal(weiz_pbar_qbar, groupColour = TRUE, groupFill = TRUE)
  
  TC_pbar_qbar <- ggplot(both_gea_merge[!(is.na(both_gea_merge$overlapping)),],aes(mean_pbar_qbar,TC_score,colour=overlapping))+
    geom_point()+
    geom_hline(yintercept=0,colour="red2")
  TC_pbar_qbar_marginal <- ggMarginal(TC_pbar_qbar, groupColour = TRUE, groupFill = TRUE)
  
  # Plot all together
  gene_summary_effect <- ggplot(both_gea_merge,aes(mean_pbar_qbar,mean_corr,colour=overlapping))+
    geom_point()
  gene_summary_effect_marginal <- ggMarginal(gene_summary_effect, groupColour = TRUE, groupFill = TRUE)
  
  # Plot all together and with SNP count
  gene_summary_effect2 <- ggplot(both_gea_merge,aes(log10(snp_count.x),mean_corr,colour=overlapping))+
    geom_point()
  gene_summary_effect_marginal2 <- ggMarginal(gene_summary_effect2, groupColour = TRUE, groupFill = TRUE)
  
  # Return figures and stats
  return(list(weiZ_TC,gene_summary_effect_marginal,gene_summary_effect_marginal2,overlap_out))
})

# Get all the combined plots and merge
pdf("figs/test_weiz_TC_by_climate.pdf",height=30,width=40)
plot_grid(plotlist = lapply(climate_TC_WZA,'[[',1))
plot_grid(plotlist = lapply(climate_TC_WZA,'[[',2))
plot_grid(plotlist = lapply(climate_TC_WZA,'[[',3))
dev.off()

# Now also pull the overlap statistics and compare these to climate stats
overlap_stats <- data.frame(rbindlist(lapply(climate_TC_WZA,'[[',3)))
climate_cline_res <- data.frame(read.table("outputs/GEA_res/Arabidopsis_halleri_Kubota_Individual/climate_cline.tsv",header = T))

# Calculate coefficient of variation for all cols
climate_var_variance <- matrix(ncol=2,nrow=length(climate_vars))
for(i in 1:length(climate_vars)){
  #climate_var_variance[i,1] <- FinCal::coefficient.variation(sd(climate_cline_res[,climate_vars[i]]),mean(climate_cline_res[,climate_vars[i]]))
  climate_var_variance[i,1] <- EnvStats::cv(climate_cline_res[,climate_vars[i]])
  climate_var_variance[i,2] <- overlap_stats[overlap_stats$var==climate_vars[i] & overlap_stats$Var1=="Both","Freq"] / sum(overlap_stats[overlap_stats$var==climate_vars[i] & overlap_stats$Var1 %in% c("Both","WZA"),"Freq"])
}

to_plot <- as.data.frame(climate_var_variance)
  ggplot(to_plot,aes(V1,V2))+
    geom_point()
