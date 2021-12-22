# Quantify peakiness
lib <- c("pbmclapply","data.table","ggplot2","viridis","ggridges","dplyr","ggExtra")
sapply(lib,library,character.only=T)

# Prepare GEA results ------------------------------------------------------
focal_datasets <- list.files("outputs/GEA_res")
all_datasets <- data.frame(rbindlist(pbmclapply(focal_datasets,function(dataset){
  # dataset="Eucalyptus_sideroxylon_Murray_Individual"
  # climate_var="mean_temp"
  
  # Fetch all the HMM results
  all_hmm <- list.files(paste0("outputs/GEA_res/",dataset,"/"),pattern="HMM_peak_clusters.rds")
  climate_vars <- gsub("_HMM_peak_clusters.rds","",all_hmm)
  
  if(length(all_hmm) == 0){
    return(NULL)
  } else {
    
    all_climates <- data.frame(rbindlist(lapply(climate_vars,function(climate_var){
      
      # Fetch HMM
      hmm_res <- readRDS(paste0("outputs/GEA_res/",dataset,"/",climate_var,"_HMM_peak_clusters.rds"))
      
      # head(hmm_res[[1]]$state_5)
      
      # pdf("~/Desktop/test.pdf",width=10,height=4)
      # for(i in 1:length(hmm_res)){
      #   g1 <- ggplot(hmm_res[[i]]$state_5,aes(pos,wza,colour=factor(state)))+
      #     geom_point()+
      #     ggtitle(names(hmm_res)[i])+
      #     theme(legend.position = "bottom")
      #   
      #   print(ggMarginal(g1,groupColour = TRUE, groupFill = TRUE,margins = "y"),newpage = TRUE)
      # }
      # dev.off()
      
      # Per chromosome calculate diff of top and middle cluster
      chrom_difs <- data.frame(rbindlist(lapply(hmm_res,function(tmp){
        
        tmp2 <- tmp$wind_100$`5_state`
        
        # cluster means
        cluster_means <- data.frame(tmp2 %>% group_by(state) %>% summarise(r2=mean(r2_mean)))
        max_state = cluster_means[order(-cluster_means$r2),"state"][1]
        
        med_state = cluster_means[which(cluster_means$r2 == median(cluster_means$r2)),"state"]
        
        # OUt
        data.frame(top_r2=cluster_means[order(-cluster_means$r2),"r2"][1],
                   four_r2=cluster_means[order(-cluster_means$r2),"r2"][2],
                   med_r2=cluster_means[order(-cluster_means$r2),"r2"][3],
                   top_winds=nrow(tmp2[tmp2$state==max_state,]),
                   med_winds=nrow(tmp2[tmp2$state==med_state,]))
      })))
      
      chrom_difs$climate=climate_var
      return(chrom_difs)
    })))
    
    all_climates$dataset=dataset
    return(all_climates)
  }
},mc.core=16)))

chrom_difs$r2_diff <- (chrom_difs$top_r2-chrom_difs$med_r2)/(chrom_difs$four_r2-chrom_difs$med_r2)
chrom_difs$r2_enrich <- chrom_difs$top_r2/chrom_difs$med_r2
chrom_difs$top_enrich <- chrom_difs$top_winds/chrom_difs$med_winds

ggplot(chrom_difs,aes(top_enrich,r2_diff))+geom_point()
ggplot(chrom_difs,aes(top_enrich,r2_enrich))+geom_point()

