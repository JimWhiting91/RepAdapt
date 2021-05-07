#######################################################################
# Pull climate data from World Clim for lat long and compare clines

# Set up environment
lib <- c("pheatmap","ggdendro","ggplot2","data.table","raster","sp","parallel","lostruct","Morpho","viridis","stars","rnaturalearth","rnaturalearthdata","cowplot")
lapply(lib,library,character.only=T)

################################################
# Function library

# Make a data.frame for climate data given lat long input and raster stack
make_climate_dataframe <- function(long_lat,climate_data){
  
  # Fetch points
  points <- SpatialPoints(long_lat, proj4string = climate_data@crs)
  values <- raster::extract(climate_data,points)
  
  # Pull data together
  df <- cbind.data.frame(coordinates(points),values)
  
  # Return the output
  return(df)
}
################################################
# Set up temp colours
temp_colors <- colorRampPalette(c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))

# Read in our metadata and filter
metadata <- read.csv("~/Calgary/RepAdapt/metadata/sample_species_vcf_author_map_v2_210427.csv")
metadata <- metadata[grep("vcf.gz",metadata$VCF),]

# ---------------
# For now, use three dummy clines, but eventually this list will be all the real clines...
# In this dummy, one cline is US, one is europe, one is Australia
# clines <- list(
#   cline1=data.frame(long=sort(sample(seq(-110,-100,0.1),10)),
#                     lat=sort(sample(seq(35,55,0.1),10))),
#   cline2=data.frame(long=sort(sample(seq(30,40,0.1),10)),
#                     lat=sort(sample(seq(50,62,0.1),10))),
#   cline3=data.frame(long=sort(sample(seq(125,145,0.1),10)),
#                     lat=sort(sample(seq(-30,-15,0.1),10)))
# )
# names(clines) <- c("USA","Europe","Australia")
# # Set up new dummy clines based on a common latitude of 15-40 in US, Asia, and Australia
# clines <- list(
#   cline1=data.frame(long=sort(sample(seq(-100,-90,0.1),50)),
#                     lat=sort(sample(seq(20,40,0.1),50))),
#   cline2=data.frame(long=sort(sample(seq(75,85,0.1),50)),
#                     lat=sort(sample(seq(20,40,0.1),50))),
#   cline3=data.frame(long=sort(sample(seq(135,145,0.1),50)),
#                     lat=sort(sample(seq(-40,-20,0.1),50)))
# )
# names(clines) <- c("USA","Asia","Australia")
# Remove this later
# ---------------

# # Fetch climate data from CLIMA dataset
# #clima_tiff <- list.files("data/climate_analysis/CLIMA/",full.names = T,pattern="tif")
# clima_tiff <- list.files("data/climate_analysis/wc10/",full.names = T,pattern="tif")
# # Create raster stack
# clima_stack <- stack(clima_tiff)

# Fetch climate data from worldclim
dir.create("data/worldclim")
clima_stack <- raster::getData(name = 'worldclim', var = 'bio', res = 2.5, path = "data/wordclim")

# Label the stack
names(clima_stack) <- c("mean_temp","mean_diurnal","isothermality","temp_seasonality","max_temp_warmest_month",
                        "min_temp_coldest_month","temp_range","mean_temp_wet_quarter","mean_temp_dry_quarter","mean_temp_warm_quarter",
                        "mean_temp_cold_quarter","annual_precip","precip_wet_month","precip_dry_month","precip_seasonality","precip_wet_quarter",
                        "precip_dry_quarter","precip_warm_quarter","precip_cold_quarter")
climate_vars <- c("Annual Mean Temperature",
                  "Mean Diurnal Range",
                  "Isothermality",
                  "Temperature Seasonality",
                  "Max Temperature of Warmest Month",
                  "Min Temperature of Coldest Month",
                  "Temperature Annual Range",
                  "Mean Temperature of Wettest Quarter",
                  "Mean Temperature of Driest Quarter",
                  "Mean Temperature of Warmest Quarter",
                  "Mean Temperature of Coldest Quarter",
                  "Annual Precipitation",
                  "Precipitation of Wettest Month",
                  "Precipitation of Driest Month",
                  "Precipitation Seasonality",
                  "Precipitation of Wettest Quarter",
                  "Precipitation of Driest Quarter",
                  "Precipitation of Warmest Quarter",
                  "Precipitation of Coldest Quarter")

# Plot each of these in ggplot
climate_figs <- mclapply(1:length(names(clima_stack)),function(x){
  
  # gplot(clima_stack[[x]]) + geom_tile(aes(fill = value)) +
  #   scale_fill_gradientn(name = climate_vars[[x]],
  #                        colors = temp_colors(5),
  #                        #    limits = c(1, 5),
  #                        na.value = "white")+
  #   coord_equal()+
  #  # ylim(-60,Inf)+
  #   scale_y_continuous(limit = c(-60, clima_stack[[1]]@extent[4])) + 
  #   theme_minimal()+
  #   theme(panel.grid = element_blank(),
  #         legend.position = "top")
  
  # Read in the climate data from raw to star
  # This code from https://bedatablog.netlify.app/post/download-and-illustrate-current-and-projected-climate-in-r/
  tmp_stars <- stars::read_stars(paste0("data/worldclim/wc10/bio",x,".bil"))
  if(x<12){
    tmp_stars <- tmp_stars/10
  }
  
  # Make the climate plot
  climate_plot <- ggplot() + 
    geom_stars(data = tmp_stars) +
    scale_fill_gradientn(name = climate_vars[x],
                         colors = temp_colors(5),
                         #limits = c(-7, 32),
                         na.value = "white") +
    coord_equal() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    theme_minimal() +
    theme(legend.position = "top",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())
  
  # Add all points
  #climate_plot <- 
  climate_plot +
    geom_point(data=metadata,aes(x=Long,y=Lat),colour="black",size=0.5,alpha=0.5)
  
},mc.cores=6)

# Plot AMT
climate_figs[[1]]

# Plot some others
plot_grid(plotlist = climate_figs[c(4,6,12,15)],
          ncol=2)

###########################################################################
#### Create our cline datasets #####
# For the timebeing, these will be based on species for simplicitys sake
species_ids <- unique(metadata$Taxon_name)
  
# Build a climate data frame for each cline
cline_data <- mclapply(1:length(species_ids),function(x){
  cat(paste0("Starting ",species_ids[x]))
  
  # Fetch all unique long-lat identifiers per species
  longlat_tmp <- na.omit(unique(metadata[metadata$Taxon_name == species_ids[x],c("Long","Lat")]))
  
  data.frame(species=species_ids[x],
             make_climate_dataframe(longlat_tmp,clima_stack))
},mc.cores=6)
names(cline_data) <- species_ids

# For each climate date.frame, we want to get PCs to compare covariance
cline_PCs <- lapply(1:length(cline_data),function(x){
  #print(x)
  
  # Fetch clim data
  clim_tmp <- na.omit(cline_data[[x]][,4:ncol(cline_data[[x]])])
  
  # Scale internally to avoid scale error...
  clim_tmp[,which(apply(clim_tmp, 2, var)==0)] <- 0
  for(col in which(apply(clim_tmp, 2, var)!=0)){
    clim_tmp[,col] <- scale(clim_tmp[,col])
  }
  
  # Make PCA  
  tmp_pca <- prcomp(clim_tmp,scale. = F,center = T)
  return(tmp_pca)
})

# For each climate date.frame, we want to get raw eigenvectors/eigenvalues for the raw correlation matrix...
# This are identical to taking the loading matrix of the PCA.
cline_eigens <- lapply(1:length(cline_PCs),function(x){
  cline_PCs[[x]]$rotation[,1]
})


##########################################################################################
# Examine the distance between centroids in common PC spaces
centroid_mat <- matrix(ncol=length(climate_vars),nrow=length(species_ids))

# Perform global PCA
global_data <- na.omit(data.frame(rbindlist(cline_data)))
global_data_clean <- global_data[,c(-1,-2,-3)]
global_pca <- prcomp(global_data_clean,scale=T,center = T)
pc_summary <- summary(global_pca)
eigenval_props <- round(pc_summary$importance[2,]*100,2)

# Fetch scores and return identifiers
global_scores <- data.frame(global_pca$x)
global_scores$species <- global_data$species

# Plot to visualise
global_climate_PCA <- ggplot(global_scores,aes(x=PC1,y=PC2,colour=species))+
  geom_point(size=2,alpha=0.75)+
  stat_ellipse(show.legend = F) +
  theme_minimal() +
  theme(
    title = element_text(size=20),
    axis.title = element_text(size=18),
    axis.text = element_text(size=16),
    legend.title = element_text(size=18),
    legend.text = element_text(size=16),
  )+
  labs(y=paste0("PC2 ",eigenval_props[2],"%"), 
       x = paste0("PC1 ",eigenval_props[1],"%"),
       color = "Species")
  #scale_color_brewer(palette="Dark2")

# Or 3D...
# library(plotly)
# plot_ly(x=global_scores[,1],y=global_scores[,2],z=global_scores[,3],type="scatter3d",mode="markers",color = global_scores[,"location"])

# Visualise the loading matrix as well...
global_loading <- data.frame(reshape2::melt(global_pca$rotation[,1:3]))
loading_fig <- ggplot(global_loading,aes(y=Var1,x=Var2,fill=value))+
  geom_tile()+
  scale_fill_gradient2()+
  theme_minimal()+
  theme(
    title = element_text(size=20),
    axis.title = element_blank(),
    axis.text = element_text(size=16),
    legend.title = element_text(size=18),
    legend.text = element_text(size=16)
  )+
  ggtitle("Global Climate PC")

# Visualise as a comparison, just mean temp and precipitation among each cline...
temp_precip <- data.frame(value=c(global_data$mean_temp,global_data$annual_precip),
                          var=rep(c("Mean Temp","Annual Precipitation"),each=nrow(global_data)),
                          species=rep(global_data$species),2)
ggplot(temp_precip,aes(y=species,x=value,fill=species))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~var,ncol=2,scales = "free_x")+
 # scale_fill_brewer(palette = "Dark2")+
  theme_bw()+
  theme(
    title = element_text(size=20),
    axis.title = element_blank(),
    axis.text = element_text(size=14),
    axis.text.x = element_text(size=14,angle=45,hjust=1),
    legend.position = "none",
    strip.text = element_text(size=18)
  )

# Get centroids
species <- unique(global_data$species)
for(i in 1:nrow(centroid_mat)){
  centroid_mat[i,] <- colMeans(global_scores[global_scores$species == species[i],grep("PC",colnames(global_scores))])
}

# And get distance mat and cluster
centroid_dists <- as.matrix(dist(centroid_mat))
rownames(centroid_dists) <- species
colnames(centroid_dists) <- species

# # Plot dendrogram...
# ddata <- dendro_data(as.dendrogram(hclust(dist(centroid_mat))), type = "rectangle")
# dist_dendro <- ggplot(segment(ddata)) + 
#   geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) + 
#   coord_flip() + 
#   scale_y_reverse(expand = c(0.2, 0))+
#   theme_dendro()
# 
# # Visualise
# # centroid_dists[lower.tri(dist(centroid_mat))] <- NA
# centroid_plot_data <- reshape2::melt(centroid_dists) 

# # Reorder according to cluster...
# ggplot(centroid_plot_data,aes(x=Var1,y=Var2,fill=value))+
#   geom_tile()+
#   theme_minimal()+
#   theme(axis.text.x = element_text(size=14,angle=45,hjust=1),
#         axis.text.y = element_text(size=14),
#         axis.title = element_blank())+
#   scale_fill_viridis(option="B")

centroid_dist_matrix <- pheatmap(centroid_dists)
centroid_dist_matrix
##########################################################################################
# Calculate the angles between the dominant eigenvectors and fill a matrix
eig_angle_matrices <- lapply(1:2,function(x){
  eigen_angles <- matrix(ncol=length(species_ids),nrow=length(species_ids))
  colnames(eigen_angles) <- names(species_ids)
  rownames(eigen_angles) <- names(species_ids)
  for(i in 1:nrow(eigen_angles)){
    for(j in 1:nrow(eigen_angles)){
      if(i != j){
        
        # Fetch our eigenvectors
        eig1 <- cline_PCs[[i]]$rotation[,x]
        eig2 <- cline_PCs[[j]]$rotation[,x]
        
        # Calculate the angle and save
        angle.deg <- angle.calc(eig1,eig2)*(180/pi)
        if(angle.deg > 90){
          angle.deg <- 180 - angle.deg
        }
        
        # Save to matrix
        eigen_angles[i,j] <- round(angle.deg,2)
      }
    }
  }
  
  # Set colnames n that
  colnames(eigen_angles) <- species_ids
  rownames(eigen_angles) <- species_ids
  
  # Reshape and return
  # eig_reshape <- reshape2::melt(eigen_angles)
  # eig_reshape$eig <- paste0("Eigenvector ",x)
  return(eigen_angles)
})

# Make a pretty heatmap
pheatmap(as.data.frame(eig_angle_matrices[[1]]))

# # Visualise the matrix...
# ggplot(eig_angle_matrices,aes(x=Var1,y=Var2,fill=value))+
#   geom_tile()+
#   facet_wrap(~eig,ncol=1,strip.position = "right")+
#   scale_fill_viridis(option="D")+
#   theme_minimal()+
#   theme(
#     title = element_text(size=20),
#     axis.title = element_blank(),
#     axis.text = element_text(size=16),
#     legend.title = element_text(size=18),
#     legend.text = element_text(size=16),
#     strip.text = element_text(size=18)
#   )+
#   labs(fill="Angle")

# # Viusalise loadings of each...
# locations <- c("USA","Asia","Australia")
# cline_loading <- NULL
# for(i in 1:length(locations)){
#   tmp <- data.frame(reshape2::melt(cline_PCs[[i]]$rotation[,1:2]))
#   tmp$location <- locations[i]
#   cline_loading <- rbind(cline_loading,tmp)
# }
# local_loading_fig <- ggplot(cline_loading,aes(y=Var1,x=Var2,fill=value))+
#   geom_tile()+
#   scale_fill_gradient2()+
#   theme_minimal()+
#   theme(
#     title = element_text(size=20),
#     axis.title = element_blank(),
#     axis.text = element_text(size=16),
#     legend.title = element_text(size=18),
#     legend.text = element_text(size=16),
#     strip.text = element_text(size=18)
#   )+
#   facet_wrap(~location,ncol=3)

##########################################################################################
# # MDS analysis based on the covariance matrices of each cline...
# 
# # We also want to look at a principal co-ordinate analysis of the loading matrices
# merged_clines <- as.matrix(rbindlist(cline_data))
# 
# # Remove the lat and long...
# merged_clines <- merged_clines[,c(-1,-2,-3)]
# class(merged_clines) <- "numeric"
# merged_clines_stan <- apply(merged_clines,2,scale)
# 
# # Use lostruct functions to perform PCA and MDS
# eigen_clines <- eigen_windows(merged_clines_stan,k=10,win=nrow(clines[[1]]))
# pcdist <- pc_dist(eigen_clines,npc=10)
# 
# # And MDS
# fit2d <- cmdscale(pcdist, eig=TRUE, k=2)
# plot( fit2d$points, xlab="Coordinate 1", ylab="Coordinate 2", col=rainbow(1.2*nrow(pcdist)) )
# 
# #-------------------------------------------------------------------------------------------------
# # WORLDCLIM DATA EXAMPLE
# # Set environmental variables
# clim <- getData("worldclim",var="bio",res=10,path = "data/climate_analysis/")
# 
# r <- r[[c(1,12)]]
# names(r) <- c("Temp","Prec")
# 
# lats <- c(9.093028 , 9.396111, 9.161417)
# lons <- c(-11.7235, -11.72975, -11.709417) 
# 
# coords <- data.frame(x=lons,y=lats)
# 
# points <- SpatialPoints(coords, proj4string = r@crs)
# 
# values <- extract(r,points)
# 
# df <- cbind.data.frame(coordinates(points),values)
# 
# df
# x        y Temp Prec
# 1 -11.72350 9.093028  257 2752
# 2 -11.72975 9.396111  257 2377
# 3 -11.70942 9.161417  257 2752