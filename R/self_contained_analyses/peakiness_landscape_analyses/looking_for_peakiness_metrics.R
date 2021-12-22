# Quantify peakiness
lib <- c("e1071","data.table","ggplot2","viridis","ggridges","dplyr","ggExtra")
sapply(lib,library,character.only=T)

# Prepare GEA results ------------------------------------------------------
focal_datasets <- list.files("outputs/GEA_res")


# Functions ---------------------------------------------------------------
plot_wza_manhattan <- function(wza_res_dir,climate_var){
  
  # Fetch the GEA res
  wza_res <- readRDS(paste0(wza_res_dir,"/",list.files(wza_res_dir,pattern=paste0(climate_var,"_WZA_TC_allgenes.rds"))))
  
  # Transform for plotting...
  wza_res$chr <- sapply(strsplit(wza_res$gene_id,":",),'[[',1)
  pos <- sapply(strsplit(wza_res$gene_id,":",),'[[',2)
  wza_res$start <- as.integer(sapply(strsplit(pos,"-",),'[[',1))
  wza_res$end <- as.integer(sapply(strsplit(pos,"-",),'[[',2))
  wza_res$mid <- rowMeans(wza_res[,c("start","end")])
  
  # Add gene index
  wza_res$gene_index <- 1:nrow(wza_res)
  
  # Plot
  ggplot(wza_res[wza_res$weiZ>0,],aes(gene_index,weiZ))+
    geom_line()+
    theme()+
    stat_density2d_filled(alpha=0.5,show.legend = F)+
    theme(axis.text=element_text(size=16),
          axis.title = element_text(size=18),
          title = element_text(size=20),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing.x=unit(0.1, "lines"),
          strip.text = element_text(angle=90,hjust=1,size=14))+
    geom_hline(yintercept = quantile(wza_res$weiZ,0.99),colour="red2")+
    geom_hline(yintercept = quantile(wza_res$weiZ[wza_res$weiZ > 0],0.5),colour="blue1")+
    labs(x="Chr/Position",y="WZA")
}

# ---------------------------------------------------------------

# Remove any focal_datasets without full wza outputs...
focal_climate="mean_temp"
focal_datasets <- na.omit(sapply(focal_datasets,function(x){
  if(file.exists(paste0("outputs/GEA_res/",x,"/",focal_climate,"_WZA_TC_allgenes.rds"))){
    return(x)
  } else {
    return(NA)
  }
}))

# Read in all of our GEA results
all_gea_res <- pbmclapply(focal_datasets,function(dataset){
  
  # Get all wza
  wza_res_tmp <- list.files(paste0("outputs/GEA_res/",dataset,"/"),pattern="_WZA_TC_allgenes.rds")
  
  all_wza <- data.frame(rbindlist(lapply(wza_res_tmp,function(file) readRDS(paste0("outputs/GEA_res/",dataset,"/",file)))))
  all_wza$dataset=dataset
  return(all_wza)
},mc.cores=n_cores)
names(all_gea_res) <- focal_datasets

# What are all combos of dataset x climate
climate_vars <- unique(all_gea_res[[1]]$climate_var)

# First loop over datasets, and quantify wza distance vs quantiles from median to max
dataset_relationships <- data.frame(rbindlist(pbmclapply(focal_datasets,function(dataset){
  
  # Subset
  tmp <- all_gea_res[[dataset]]
  
  # Now loop over climate vars
  clim_peakiness <- data.frame(rbindlist(lapply(climate_vars,function(var){
    
    # Subset again
    tmp2 <- tmp[tmp$climate_var==var,]
    
    # Take 10 quantiles between median and max
    quants <- seq(0.5,1,by=0.01)
    
    # Get wza of each...
    quant_wza <- quantile(tmp2$weiZ,quants)
    
    # Now calculate the distance between points...
    out <- data.frame(quants=quants[2:length(quants)],
                      wza=quant_wza[2:length(quants)],
                      climate_var=var)
    
    for(i in 1:nrow(out)){
      out$dist[i] <- out$wza[i]-quant_wza[i]
    }
    return(out)
  })))
  
  # Visualise by climate variable...
  clim_peakiness$dataset=dataset
  
  
  return(clim_peakiness)
},mc.cores=4)))

# Visualise the first 10...
ggplot(dataset_relationships[dataset_relationships$dataset %in% focal_datasets[1:10],],aes(x=quants,y=dist,colour=climate_var))+
  geom_line()+
  facet_wrap(~dataset)+
  scale_y_continuous(trans="log")

# Also fetch skew of the right hand side of wzas
all_gea_res_bind <- data.frame(rbindlist(all_gea_res))
skew_res <- data.frame(all_gea_res_bind %>% group_by(climate_var,dataset) %>% summarise(skew=skewness(weiZ),
                                                                                        kurtosis=kurtosis(weiZ)))
skew_res$climate_dataset <- paste0(skew_res$climate_var,"-",skew_res$dataset)

ggplot(skew_res,aes(x=skew,y=kurtosis))+geom_point()

# Save skew res
saveRDS(skew_res,"outputs/all_GEA_skew_results_tmp.rds")

# Order
# Low
tail(skew_res[order(-skew_res$skew),],40)
plot_wza_manhattan("outputs/GEA_res/Arabidopsis_lyrata_Savolainen_Individual","temp_range")

# Low
head(skew_res[order(-skew_res$skew),],10)
plot_wza_manhattan("outputs/GEA_res/Eucalyptus_sideroxylon_Murray_Individual","mean_temp")
plot_wza_manhattan("outputs/GEA_res/Helianthus_argophyllus_Todesco_Individual","mean_temp")

# Save the skew results somewhere...
saveRDS(skew_res,"outputs/wza_skewness.rds")


# Quantify skew decay -----------------------------------------------------

all_climate_decay_zero <- data.frame(rbindlist(lapply(climate_vars,function(climate_var){
  all_decay_zero <- sapply(focal_datasets,function(dataset){
    eucalypt_meantemp <- all_gea_res_bind[all_gea_res_bind$dataset==dataset &
                                            all_gea_res_bind$climate_var==climate_var,]
    
    eucalypt_meantemp <- eucalypt_meantemp[order(-eucalypt_meantemp$weiZ),]
    
    skew_decay <- data.frame(rbindlist(pbmclapply(0:1000,function(i){
      
      out=data.frame(removed=i,
                     skew=skewness(eucalypt_meantemp$weiZ[(i+1):nrow(eucalypt_meantemp)]))
      
    },mc.cores=4)))
    
    skew_decay$log_skew <- log(skew_decay$skew)
    nrow(na.omit(skew_decay))
  })
  
  to_plot <- data.frame(skew_zero = all_decay_zero,
                        dataset = focal_datasets,
                        climate_var = climate_var)
  return(to_plot)
})))

all_climate_decay_zero$climate_dataset <- paste0(all_climate_decay_zero$climate_var,"-",all_climate_decay_zero$dataset)
to_plot <- merge(all_climate_decay_zero,skew_res[,c("skew","kurtosis","climate_dataset")])
ggplot(to_plot,aes(x=skew,y=skew_zero))+
  geom_point()

head(all_climate_decay_zero[order(-all_climate_decay_zero$skew_zero),],20)

ggplot(all_gea_res_bind[all_gea_res_bind$dataset %in% c("Helianthus_petiolaris_Todesco_Individual","Eucalyptus_albens_Murray_Individual"),],aes(x=weiZ,y=dataset))+
  geom_density_ridges()+
  xlim(10,100)

## Build a null based on rnorm with skew? Predict the analytical null?
## Kurtosis

