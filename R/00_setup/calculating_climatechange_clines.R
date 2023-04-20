library(raster)
library(sp)
library(pbmcapply)

source("R/repadapt_functions.R")

# Fetch all the climate rasters...
# Precipitation...
raster_dir = "data/worldclim/wc2.1_2.5m_prec_1960-1969/"
prec_rasters = paste0(raster_dir,"/",list.files(raster_dir,pattern = "tif"))
raster_dir = "data/worldclim/wc2.1_2.5m_prec_2010-2018/"
prec_rasters2 = paste0(raster_dir,"/",list.files(raster_dir,pattern = "tif"))
prec_rasters = c(prec_rasters,prec_rasters2)

# Max Temperature...
raster_dir = "data/worldclim/wc2.1_2.5m_tmax_1960-1969/"
tmax_rasters = paste0(raster_dir,"/",list.files(raster_dir,pattern = "tif"))
raster_dir = "data/worldclim/wc2.1_2.5m_tmax_2010-2018/"
tmax_rasters2 = paste0(raster_dir,"/",list.files(raster_dir,pattern = "tif"))
tmax_rasters = c(tmax_rasters,tmax_rasters2)

# Find all sampling files and make climate change files -------------------
sampling_data_list = list.files("data/VCFs/",pattern = "sampling_data.csv", recursive = TRUE,full.names = T)
sampling_data_list = grep('old',sampling_data_list,invert = T,value = T)

# Loop over these and make a climate_change_env.txt file
for(sampling_input in sampling_data_list){
  print(paste0(">> STARTING ",which(sampling_data_list == sampling_input)," of ",length(sampling_data_list)))
  
  if(file.exists(paste0(dirname(sampling_input),"/climate_change_env.txt"))){
    print("Climate change data already exists, skipping...")
  } else {
    
    clim_change_out = calculate_climchange_from_longlat(sampling_input)
    write.table(clim_change_out,
                paste0(dirname(sampling_input),"/climate_change_env.txt"),
                quote = F,row.names = F,sep = "\t")
  }
}

# How do climate change rasters associate with other clim vars... ---------
# For every dataset, run through and get the correlation coefficient of climate change vars vs regular climate vars...
data_dirs = list.files("outputs/GEA_res/",pattern = "220927")
data_dirs = grep(".rds",data_dirs,invert = T,value = T)
clim_change_env_corrs = rbindlist(lapply(data_dirs,function(res){
  print(res)
  clim_cline = read.table(paste0("outputs/GEA_res/",res,"/climate_cline.tsv"),header=T)
  
  corr_tmp = melt(cor(clim_cline[,3:ncol(clim_cline)],method = "spearman"))
  tmax_corrs = corr_tmp[corr_tmp$Var1 == "tmax_clim_change" & !(corr_tmp$Var2 %in% c("tmax_clim_change","prec_clim_change")),]
  prec_corrs = corr_tmp[corr_tmp$Var1 == "prec_clim_change" & !(corr_tmp$Var2 %in% c("tmax_clim_change","prec_clim_change")),]
  
  out = data.table(rbind(tmax_corrs,prec_corrs),res = res)
}))

# Fetch the bioclim vars
bioclim = colnames(read.table(paste0("outputs/GEA_res/",data_dirs[1],"/climate_cline.tsv"),header=T))[3:21]

# Plot these as correlation distributions...
clim_change_env_corrs$Var2 = stringr::str_to_title(gsub("_"," ",clim_change_env_corrs$Var2))
clim_change_env_corrs$clim_change_F = stringr::str_to_title(gsub("_"," ",clim_change_env_corrs$Var1))
clim_change_env_corrs$clim_F = factor(clim_change_env_corrs$Var2,levels = stringr::str_to_title(gsub("_"," ",bioclim)))
clim_change_clim_env_corrs = ggplot(clim_change_env_corrs, aes(y = clim_F,fill = clim_change_F,x = value)) +
  geom_boxplot() +
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text = element_text(size = 12),
        axis.title.x = element_text(size = 14),
        legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 11)) +
  labs(x = expression(Spearmans~rho),fill = "") +
  geom_vline(xintercept = 0,linetype = "dotted") +
  scale_fill_brewer(palette = "Set1")

pdf("figs/FigureSX_corr_between_climchange_and_climenv.pdf",height=8,width = 6)
clim_change_clim_env_corrs
dev.off()

