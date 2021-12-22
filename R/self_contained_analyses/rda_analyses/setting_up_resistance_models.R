### Build a raster object for circuitscape modelling...
lib <- c("stars","ggplot2","raster","data.table","dplyr")
sapply(lib,library,character.only=T)

# Test with Amaranthus dataset...
dataset = "Amaranthus_tuberculatus_Wright_Individual"
cline_input <- read.table(paste0("outputs/GEA_res/",dataset,"/climate_cline.tsv"),header=T)

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

# Set climate labs
climate_labs <- c("mean_temp",
                  "mean_diurnal",
                  "isothermality",
                  "temp_seasonality",
                  "max_temp_warmest_month",
                  "min_temp_coldest_month",
                  "temp_range",
                  "mean_temp_wet_quarter",
                  "mean_temp_dry_quarter",
                  "mean_temp_warm_quarter",
                  "mean_temp_cold_quarter",
                  "annual_precip",
                  "precip_wet_month",
                  "precip_dry_month",
                  "precip_seasonality",
                  "precip_wet_quarter",
                  "precip_dry_quarter",
                  "precip_warm_quarter",
                  "precip_cold_quarter")
names(climate_vars) <- climate_labs

# Set up temp colours
temp_colors <- colorRampPalette(c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))

# Read in desired raster...
climate_input <- 2 # mean_diurnal
tmp_stars <- read_stars(paste0("data/worldclim/wc2-5/bio",climate_input,".bil"))

# Set the limits
long_range <- max(cline_input$Long) - min(cline_input$Long)
long_min <- min(cline_input$Long) - 0.2*long_range
long_max <- max(cline_input$Long) + 0.2*long_range

lat_range <- max(cline_input$Lat) - min(cline_input$Lat)
lat_min <- min(cline_input$Lat) - 0.2*lat_range
lat_max <- max(cline_input$Lat) + 0.2*lat_range

# Filter the stars
tmp_stars_part <- stars:::filter.stars(tmp_stars, x > long_min, x < long_max,y > lat_min,y < lat_max)

# Plot to check
ggplot() + 
  geom_stars(data = tmp_stars_part) +
  scale_fill_gradientn(name = climate_vars[climate_input],
                       colors = temp_colors(5),
                       #limits = c(-7, 32),
                       na.value = "white") +
  coord_equal() +
  # scale_x_discrete(expand = c(0, 0)) +
  # scale_y_discrete(expand = c(0, 0)) +
  theme(legend.position = "top",
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank())+
  scale_x_continuous(limits = c(long_min,long_max))+
  scale_y_continuous(limits = c(lat_min,lat_max))+
  geom_point(data=cline_input,aes(x=Long,y=Lat),colour="black",size=2,alpha=1)

# Make resistance map...
# Header
r_header <- matrix(nrow=6,ncol=2)
r_header[,1] <- c("ncols","nrows","xllcorner","yllcorner","cellsize","NODATA_value")
r_header[,2] <- c(nrow(tmp_stars_part[[1]]),ncol(tmp_stars_part[[1]]),1,1,1,-9999)

# Body
r_body <- t(as.matrix(tmp_stars_part[[1]]))
r_body[is.na(r_body)] <- -9999

# # ugly write and merge
# write.table(r_header,"outputs/header_tmp",quote = F,sep = "\t",row.names = F,col.names = F)
# write.table(r_body,"outputs/body_tmp",quote = F,sep = "\t",row.names = F,col.names = F)
# system("cat outputs/header_tmp outputs/body_tmp > outputs/kubota_Ahalleri_mean_diurnal_raster2.tif")
# system("rm -f outputs/header_tmp outputs/body_tmp")
# 
# # And write locations file...
# long_seq <- seq(long_min,long_max,(long_max-long_min)/(ncol(r_body)-1))
# lat_seq <- rev(seq(lat_min,lat_max,(lat_max-lat_min)/(nrow(r_body)-1)))
# locations_mat <- matrix(nrow=nrow(climate_cline),ncol = 3)
# locations_mat[,1] <- 1:nrow(locations_mat)
# for(i in 1:nrow(locations_mat)){
#   locations_mat[i,2] <- length(long_seq[long_seq < climate_cline$Long[i]])+1.5
#   locations_mat[i,3] <- length(lat_seq[lat_seq < climate_cline$Lat[i]])+1.5
# }
# 
# write.table(locations_mat,
#             "outputs/kubota_Ahalleri_resistance_positions.txt",
#             sep="\t",col.names = F,quote = F,row.names = F)
# 
# # Plotting test
# points_plot <- data.frame(locations_mat)
# 
# # Plot body to check locations match up with lat long
# r_body_tmp <- r_body
# body_melt <- reshape2::melt(r_body_tmp)
# body_melt$value[body_melt$value == -9999] <- NA
# ggplot(body_melt,aes(Var2,rev(Var1)))+
#   geom_tile(aes(fill=value))+
#   geom_point(data=points_plot,aes(x=X2,y=X3))

# # Write
# write_stars(tmp_stars_part,'outputs/kubota_Ahalleri_mean_diurnal_raster.tif',layer="bio2.bil")
# 
# test <- st_as_sf(tmp_stars_part,merge=T)
# 
# # Merge here merges identical values to a single cell
# test <- st_as_sf(tmp_stars_part, as_points = FALSE, merge = FALSE)
# plot(test)
# st_write(test,"outputs/kubota_Ahalleri_mean_diurnal_raster.tif")
# 
# # Rasterise
# r <- 
# values(r) <- sample(x = 1:10,size = ncell(r),replace = T)
# 
# writeRaster(r,'test.tif',options=c('TFW=YES'))
# tmp_raster
# 
# # Save raster object...
# # Convert points to sp (assumes that the sf object is called example_points)
# example_points <- as(st_as_sf(tmp_stars_part,merge=T), "Spatial")
# 
# # Generate empty raster layer and rasterize points
# example_raster <- raster(crs = crs(example_points), vals = 0, resolution = c(0.5, 0.5), ext = extent(c(long_min, long_max, lat_min, lat_max))) %>%
#   rasterize(example_points, .)
# 
# writeRaster(example_raster,'outputs/kubota_Ahalleri_mean_diurnal_raster.tif',overwrite=T)
# 


# Testing Gdistance approach... -------------------------------------------
library(raster)
library(gdistance)

# Build raster from the above...
# Convert points to sp (assumes that the sf object is called example_points)
test_points <- as(st_as_sf(tmp_stars_part,merge=F), "Spatial")

# Generate empty raster layer and rasterize points
test_raster <- raster(ncol=ncol(r_body), nrow=nrow(r_body),ext = extent(c(long_min, long_max, lat_min, lat_max)))
r.polys <- rasterize(test_points, test_raster, field = test_points@data[,1], fun = "mean", 
                     update = TRUE, updateValue = "NA")
plot(r.polys)

# Convert our sampling locations to a raster spatial points object and retain projection
sampling_pos <- SpatialPointsDataFrame(coords = cline_input[,c("Long","Lat")],data=cline_input,
                               proj4string = crs(r.polys))
points(sampling_pos, pch=3)

# Calculate conductance values...
layer_cost <- 1/r.polys$layer

# Create a transition layer...
tr.cost1 <- gdistance::transition(layer_cost, transitionFunction=mean, directions=8) 
raster::plot(raster::raster(tr.cost1))
points(sampling_pos, pch=3)

# Test this raster...
# writeRaster(raster::raster(tr.cost1),'outputs/kubota_Ahalleri_mean_diurnal_raster.asc',overwrite=T,format="ascii")


# Correct for geometric distortion - type "r" corrects for lat-long distortion for random walks...
tr.cost1 <- gdistance::geoCorrection(tr.cost1,type = "r",multpl=FALSE)

# Calculate physical distance
phys.dist <- 

# Calculate least-cost distance
cost1.dist <- gdistance::costDistance(tr.cost1,sampling_pos[,1:2])

# Calculate commute time (analogous to circuitscape...)
comm1.dist <- gdistance::commuteDistance(x = tr.cost1, coords = sampling_pos[,1:2])

# Compare them
to_plot <- data.frame(LCP=as.numeric(cost1.dist),
                      COMM=as.numeric(comm1.dist))
plot(to_plot$COMM,test$COMM)
ggplot(to_plot,aes(LCP,COMM))+
  geom_point()

# # Calculate pairwise 'passage'...
# resistance_mat <- matrix(ncol=length(sampling_pos),nrow=length(sampling_pos))
# for(i in 1:length(sampling_pos)){
#   print(paste0("Starting ",i," of ",length(sampling_pos)))
#   for(j in 1:sampling_pos){
#     if(i > j){
#       res_dist <- passage(tr.cost1,sampling_pos[i,1:2],sampling_pos[j,1:2])
#     }
#   }
# }

#


# -------------------------------------------------------------------------

err.cost <- (1/RasterMaps$err27)
ffp.cost <- (RasterMaps$ffp/5)
gsp.cost <- (RasterMaps$gsp-196)/15
cti.cost <- RasterMaps$cti/5
cost1 <- (gsp.cost + cti.cost + err.cost + ffp.cost)
plot(cost1)
tr.cost1 <- gdistance::transition(cost1, transitionFunction=mean, directions=8) 
tr.cost1 <- gdistance::geoCorrection(tr.cost1,type = "c",multpl=FALSE)
comm1.dist <- gdistance::commuteDistance(x = tr.cost1, coords = sites)

