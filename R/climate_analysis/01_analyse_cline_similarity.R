#######################################################################
# Pull climate data from World Clim for lat long and compare clines

# Set up environment
lib <- c("ggplot2","data.table","raster","sp","parallel","lostruct","Morpho","viridis")
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
# Set up new dummy clines based on a common latitude of 15-40 in US, Asia, and Australia
clines <- list(
  cline1=data.frame(long=sort(sample(seq(-100,-90,0.1),50)),
                    lat=sort(sample(seq(20,40,0.1),50))),
  cline2=data.frame(long=sort(sample(seq(75,85,0.1),50)),
                    lat=sort(sample(seq(20,40,0.1),50))),
  cline3=data.frame(long=sort(sample(seq(135,145,0.1),50)),
                    lat=sort(sample(seq(-40,-20,0.1),50)))
)
names(clines) <- c("USA","Asia","Australia")
# Remove this later
# ---------------

# Fetch climate data from CLIMA dataset
clima_tiff <- list.files("data/climate_analysis/CLIMA/",full.names = T,pattern="tif")

# Create raster stack
clima_stack <- stack(clima_tiff)

names(clima_stack) <- c("mean_temp","mean_diurnal","isothermality","temp_seasonality","max_temp_warmest_month",
                        "min_temp_coldest_month","temp_range","mean_temp_wet_quarter","mean_temp_dry_quarter","mean_temp_warm_quarter",
                        "mean_temp_cold_quarter","annual_precip","precip_wet_month","precip_dry_month","precip_seasonality","precip_wet_quarter",
                        "precip_dry_quarter","precip_warm_quarter","precip_cold_quarter")

# Visualise some of these...
plot(clima_stack[[1]]) + title("Mean Annual Temp")
plot(clima_stack[[4]]) + title("Mean Temp Seasonality")
plot(clima_stack[[6]]) + title("Min Coldest Temp")

plot(clima_stack[[11]]) + title("Mean Precipitation")
plot(clima_stack[[14]]) + title("Mean Precipitation Seasonality")
plot(clima_stack[[13]]) + title("Min Precipitation Driest Month")

# Build a climate data frame for each cline
cline_data <- mclapply(1:length(clines),function(x){
  data.frame(location=names(clines)[x],make_climate_dataframe(clines[[x]],clima_stack))
},mc.cores=6)

# For each climate date.frame, we want to get PCs to compare covariance
cline_PCs <- lapply(cline_data,function(x){
  tmp_pca <- prcomp(na.omit(x[,3:ncol(x)]),scale. = T,center = T)
  return(tmp_pca)
})

# For each climate date.frame, we want to get raw eigenvectors/eigenvalues for the raw correlation matrix...
# This are identical to taking the loading matrix of the PCA.
cline_eigens <- lapply(cline_data,function(x){
  tmp_cor <- cor(na.omit(x[,3:ncol(x)]))
  tmp_eigen <- eigen(tmp_cor)
  return(tmp_eigen)
})


##########################################################################################
# Examine the distance between centroids in common PC spaces
centroid_matrix <- matrix(ncol=length(clines),nrow=length(clines))

# Perform global PCA
global_data <- na.omit(data.frame(rbindlist(cline_data)))
global_data_clean <- global_data[,c(-1,-2,-3)]
global_pca <- prcomp(global_data_clean,scale=T,center = T)
pc_summary <- summary(global_pca)
eigenval_props <- round(pc_summary$importance[2,]*100,2)

# Fetch scores and return identifiers
global_scores <- data.frame(global_pca$x)
global_scores$location <- global_data$location

# Plot to visualise
global_climate_PCA <- ggplot(global_scores,aes(x=PC1,y=PC2,colour=location))+
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
       color = "Location")+
  scale_color_brewer(palette="Dark2")+
  ggtitle("Climate over lat 20-40 degrees")

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
                          location=rep(global_data$location),2)
ggplot(temp_precip,aes(x=location,y=value,fill=location))+
  geom_boxplot(outlier.shape = NA)+
  facet_wrap(~var,ncol=2,scales = "free_y")+
  scale_fill_brewer(palette = "Dark2")+
  theme_minimal()+
  theme(
    title = element_text(size=20),
    axis.title = element_blank(),
    axis.text = element_text(size=14),
    axis.text.x = element_text(size=14,angle=45,hjust=1),
    legend.position = "none",
    strip.text = element_text(size=18)
  )

# Get centroids
locations <- unique( global_data$location)
for(i in 1:nrow(centroid_mat)){
  centroid_mat[i,] <- colMeans(global_scores[global_scores$location == locations[i],grep("PC",colnames(global_scores))])
}

# And get distance mat
centoid_dists <- dist(centroid_mat)

##########################################################################################
# Calculate the angles between the dominant eigenvectors and fill a matrix
eig_angle_matrices <- data.frame(rbindlist(lapply(1:2,function(x){
  eigen_angles <- matrix(ncol=length(clines),nrow=length(clines))
  colnames(eigen_angles) <- names(clines)
  rownames(eigen_angles) <- names(clines)
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
  
  # Reshape and return
  eig_reshape <- reshape2::melt(eigen_angles)
  eig_reshape$eig <- paste0("Eigenvector ",x)
  return(eig_reshape)
})))

# Visualise the matrix...
ggplot(eig_angle_matrices,aes(x=Var1,y=Var2,fill=value))+
  geom_tile()+
  facet_wrap(~eig,ncol=1,strip.position = "right")+
  scale_fill_viridis(option="D")+
  theme_minimal()+
  theme(
    title = element_text(size=20),
    axis.title = element_blank(),
    axis.text = element_text(size=16),
    legend.title = element_text(size=18),
    legend.text = element_text(size=16),
    strip.text = element_text(size=18)
  )+
  labs(fill="Angle")

# Viusalise loadings of each...
locations <- c("USA","Asia","Australia")
cline_loading <- NULL
for(i in 1:length(locations)){
  tmp <- data.frame(reshape2::melt(cline_PCs[[i]]$rotation[,1:2]))
  tmp$location <- locations[i]
  cline_loading <- rbind(cline_loading,tmp)
}
local_loading_fig <- ggplot(cline_loading,aes(y=Var1,x=Var2,fill=value))+
  geom_tile()+
  scale_fill_gradient2()+
  theme_minimal()+
  theme(
    title = element_text(size=20),
    axis.title = element_blank(),
    axis.text = element_text(size=16),
    legend.title = element_text(size=18),
    legend.text = element_text(size=16),
    strip.text = element_text(size=18)
  )+
  facet_wrap(~location,ncol=3)

##########################################################################################
# MDS analysis based on the covariance matrices of each cline...

# We also want to look at a principal co-ordinate analysis of the loading matrices
merged_clines <- as.matrix(rbindlist(cline_data))

# Remove the lat and long...
merged_clines <- merged_clines[,c(-1,-2,-3)]
class(merged_clines) <- "numeric"
merged_clines_stan <- apply(merged_clines,2,scale)

# Use lostruct functions to perform PCA and MDS
eigen_clines <- eigen_windows(merged_clines_stan,k=10,win=nrow(clines[[1]]))
pcdist <- pc_dist(eigen_clines,npc=10)

# And MDS
fit2d <- cmdscale(pcdist, eig=TRUE, k=2)
plot( fit2d$points, xlab="Coordinate 1", ylab="Coordinate 2", col=rainbow(1.2*nrow(pcdist)) )

#-------------------------------------------------------------------------------------------------
# WORLDCLIM DATA EXAMPLE
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