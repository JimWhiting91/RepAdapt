# Pulling gbif occurrences to get at global climate niches
lib = c('data.table',
        'dplyr',
        'ggplot2',
        'rgbif',
        'raster',
        'sf',
        'sp',
        'dynRB',
        'geosphere',
        'rnaturalearth')
sapply(lib,library,character.only = T)

# Make a data.frame for climate data given lat long input and raster stack
make_climate_dataframe <- function(long_lat,climate_data){
  
  # Fetch points
  points <- SpatialPoints(long_lat, proj4string = climate_data@crs)
  values <- terra::extract(climate_data,points)
  
  # Pull data together
  df <- cbind.data.frame(coordinates(points),values)
  
  # Return the output
  return(df)
}

# Set some variables and fetch the GEA results
run_name = "230321"
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"
pvals_file = paste0("outputs/GEA_res/run",run_name,"_",output_name,"_WZA_OG_PerGene_pvals.rds")
OG_pergene_pvals = readRDS(pvals_file)

# Loop through the species data and fetch all occurrences with lat lon
all_species = unique(OG_pergene_pvals$mean_temp$species)
# Remove the hybrid conifer and search for the parental species
# Will combine these at the end
all_species = all_species[all_species != 'Picea glaucaxengelmannii']
all_species = c(all_species,'Picea glauca','Picea engelmannii')

# Get the occ counts for all species
all_species_occ_counts = sapply(all_species,function(species) occ_count(scientificName = species))

# Get all the taxonkeys
species_taxonKeys = sapply(all_species,function(x){
  tmp = name_backbone(name = x)
  tmp$usageKey
})

# Use occ_download to pull all data for all species
if(!file.exists(paste0('outputs/',run_name,'_species_gbif_occs.rds'))){
  occ_download(
    type="and",
    pred("hasGeospatialIssue", FALSE),
    pred("hasCoordinate", TRUE),
    pred("occurrenceStatus","PRESENT"), 
    pred_or(  
      pred_lt("coordinateUncertaintyInMeters",10000),
      pred_isnull("coordinateUncertaintyInMeters")
    ),
    pred_in("taxonKey", species_taxonKeys),
    format = "SIMPLE_CSV"
  )
  
  # <<gbif download>>
  # Created: 2024-05-03T09:30:29.394+00:00
  # Citation Info:  
  #   Please always cite the download DOI when using this data.
  # https://www.gbif.org/citation-guidelines
  # DOI: 10.15468/dl.mcbger
  # Citation:
  #   GBIF Occurrence Download https://doi.org/10.15468/dl.mcbger Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-05-03
  
  d <- occ_download_get('0013915-240425142415019') %>%
    occ_download_import() |>
    as.data.table()
  
  # Only keep useful co-ords
  d_coord = d[,.(species,decimalLatitude,decimalLongitude)]
  # How many observations?
  table(d_coord$species)
  
  # Save these
  saveRDS(d_coord,
          paste0('outputs/',run_name,'_species_gbif_occs.rds'))
}
d_coord = readRDS(paste0('outputs/',run_name,'_species_gbif_occs.rds'))


# Get a measure of average distance between occs -----------------------
# Will do this only for within-continent comparisons, to exclude large
# cross-ocean distances

# Get the continent for every occurrence...
world = ne_countries(scale = "medium", returnclass = "sf")
# Check for invalid geometries and attempt to repair them
world = st_make_valid(world)
# Simplify geometries to avoid complex operations errors
world = st_simplify(world, preserveTopology = TRUE)

# Transform and spatial join
all_occs = st_as_sf(d_coord, coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)
locations_with_continents = st_join(all_occs, world, join = st_within)

# Save this as a data table
all_occs_continent = as.data.table(locations_with_continents) |>
  dplyr::select(species,continent)
all_occs_continent = cbind(all_occs_continent,st_coordinates(locations_with_continents)) |>
  dplyr::rename(lon = X,
                lat = Y)
all_occs_continent = na.omit(all_occs_continent)

# Calculate some distances...
set.seed(1234)
distance_estimates = lapply(unique(all_occs_continent$species),function(sp){
  print(paste0('>>> Starting ',sp))
  
  # Count up number of occs
  sp_occs = all_occs_continent[species == sp,]
  
  # If >1,0000, do 100 draws of 1,000
  if(nrow(sp_occs) > 1000){
    # Draw 1,000
    rand_draws = pbmcapply::pbmclapply(1:100,function(iter){
      sp_occs_sub = sp_occs[sample(1:nrow(sp_occs),1000,replace = F),]
      
      # Calculate distance matrix
      distance_matrix <- distm(sp_occs_sub[,.(lon,lat)], 
                               fun = distHaversine) / 1000  # Convert meters to kilometers
      # Calculate the mean within-continent distances
      cont_dists = sapply(unique(sp_occs_sub$continent),function(c){
        mean(distance_matrix[which(sp_occs_sub$continent == c),
                             which(sp_occs_sub$continent == c)])
      })
      weighted.mean(cont_dists,
                    w = table(sp_occs_sub$continent)[names(cont_dists)]/nrow(sp_occs_sub))
    },mc.cores = 6)
    return(mean(unlist(rand_draws)))
  } else {
    
    # Just do the same for the occs we have...
    # Calculate distance matrix
    distance_matrix <- distm(sp_occs[,.(lon,lat)], 
                             fun = distHaversine) / 1000  # Convert meters to kilometers
    # Calculate the mean within-continent distances
    cont_dists = sapply(unique(sp_occs$continent),function(c){
      mean(distance_matrix[which(sp_occs$continent == c),
                           which(sp_occs$continent == c)])
    })
    out = weighted.mean(cont_dists,
                  w = table(sp_occs$continent)[names(cont_dists)]/nrow(sp_occs))
    return(out)
  }
  
  
})

# Average through the hybrid species and return that
glaucaxengelmanni_dist = distance_estimates[species %in% c('Picea glauca',
                                                           'Picea engelmannii'),avg_gbif_distance] |>
  mean()
distance_estimates = rbind(distance_estimates[!species %in% c('Picea glauca','Picea engelmannii'),],
                           data.table(species = 'Picea glaucaxengelmannii',
                                      avg_gbif_distance = glaucaxengelmanni_dist))

# Now need to calculate the sampled distances in the same way...
# Fetch all of the original climate clines
res_dir = grep(".rds",list.files("outputs/GEA_res/",run_name),invert = T,value = T)

sampled_clines = rbindlist(lapply(paste0("outputs/GEA_res/",res_dir,"/climate_cline.tsv"),function(x){
  tmp = fread(x)
  tmp$species = paste(strsplit(gsub(paste0(run_name,'_'),'',basename(dirname(x))),
                               '_')[[1]][1:2],
                      collapse = ' ')
  tmp$dataset = basename(dirname(x))
  tmp
}))

# Get the continents...
sampled_occs = st_as_sf(sampled_clines, coords = c("Long", "Lat"), crs = 4326)
sampled_locations_with_continents = st_join(sampled_occs, world, join = st_within)

# Save this as a data table
sampled_occs_continent = as.data.table(sampled_locations_with_continents) |>
  dplyr::select(species,continent,dataset)
sampled_occs_continent = cbind(sampled_occs_continent,
                               st_coordinates(sampled_locations_with_continents)) |>
  dplyr::rename(lon = X,
                lat = Y)
sampled_occs_continent = na.omit(sampled_occs_continent)

sampled_distance_estimates = lapply(unique(sampled_occs_continent$dataset),function(d){
  print(paste0('>>> Starting ',d))
  
  # Count up number of occs
  sp_occs = sampled_occs_continent[dataset == d,]
  
  # Calculate distance matrix
  distance_matrix <- distm(sp_occs[,.(lon,lat)], 
                           fun = distHaversine) / 1000  # Convert meters to kilometers
  # Calculate the mean within-continent distances
  cont_dists = sapply(unique(sp_occs$continent),function(c){
    mean(distance_matrix[which(sp_occs$continent == c),
                         which(sp_occs$continent == c)])
  })
  out = weighted.mean(cont_dists,
                      w = table(sp_occs$continent)[names(cont_dists)]/nrow(sp_occs))
  
  return(data.table(avg_sampled_dist = out,
                    species = unique(sp_occs$species),
                    dataset = d))
  
}) |>
  rbindlist()

# Avg within species
sampled_distance_estimates = sampled_distance_estimates[,.(avg_sampled_dist = mean(avg_sampled_dist)),by = species]
sampled_distance_estimates = merge(x = sampled_distance_estimates,y = distance_estimates)

# Save these
write.csv(sampled_distance_estimates,
          paste0('outputs/',output_name,'_gbif_avg_distance.csv'),
          quote = F,row.names = F)

# Build bioclim niches for each species -----------------------------------
message(">>> Assembling Climate Clines Based on WorldClim 2.5")

# Fetch the climate data
dir.create("data/climate_data",showWarnings = F)
dir.create("data/climate_data/wordclim",showWarnings = F)
clima_stack <- raster::getData(name = 'worldclim', var = 'bio', res = 2.5, path = "data/climate_data/wordclim")

# Label the stack
names(clima_stack) <- c("mean_temp","mean_diurnal","isothermality","temp_seasonality","max_temp_warmest_month",
                        "min_temp_coldest_month","temp_range","mean_temp_wet_quarter","mean_temp_dry_quarter","mean_temp_warm_quarter",
                        "mean_temp_cold_quarter","annual_precip","precip_wet_month","precip_dry_month","precip_seasonality","precip_wet_quarter",
                        "precip_dry_quarter","precip_warm_quarter","precip_cold_quarter")

# Split by species and parallelise
species_niches = pbmcapply::pbmclapply(unique(d_coord$species),function(sp){
  
  # Fetch the env data...
  species_niche <- make_climate_dataframe(d_coord[species == sp,.(decimalLongitude,decimalLatitude)],
                                          climate_data = clima_stack)
  species_niche = na.omit(species_niche)
  species_niche$species = sp
  species_niche
},mc.cores = 6)
names(species_niches) = unique(d_coord$species)

# Combine back the hybrid species...
glauca_engelmannii = rbind(species_niches$`Picea glauca`,
                           species_niches$`Picea engelmannii`)
glauca_engelmannii$species = 'Picea glaucaxengelmannii'
species_niches = species_niches[!names(species_niches) %in% c('Picea glauca','Picea engelmannii')]
species_niches$`Picea glaucaxengelmannii` = glauca_engelmannii

# The global climate niche for each species is described by a PCA done within each species across
# all GBIF data
species_niche_pca = lapply(species_niches,function(x){
  print(x$species[1])
  prcomp(x[,!colnames(x) %in% c('species','decimalLongitude','decimalLatitude')],
         scale. = T,
         center = T)
})
names(species_niche_pca) = names(species_niches)

# Compare hypervolume overlap of sampled distributions --------------------
# Fetch all of the original climate clines
res_dir = grep(".rds",list.files("outputs/GEA_res/",run_name),invert = T,value = T)

sampled_clines = rbindlist(lapply(paste0("outputs/GEA_res/",res_dir,"/climate_cline.tsv"),function(x){
  tmp = fread(x)
  tmp$species = paste(strsplit(gsub(paste0(run_name,'_'),'',basename(dirname(x))),
                               '_')[[1]][1:2],
                      collapse = ' ')
  tmp
}))

# Now run through each species and calculate the hypervolume size and overlap
# Project the sampled clines into the PC space defined by the global dataset...
species_hypervolume_overlap = pbmcapply::pbmclapply(unique(sampled_clines$species),function(sp){
  
  # Get the projected PC scores
  trans_cline = data.frame(sampled_clines)[sampled_clines$species == sp,rownames(species_niche_pca[[sp]]$rotation)]
  # Process these
  for(var in names(species_niche_pca[[sp]]$scale)){
    trans_cline[,var] = scale(trans_cline[,var],
                              center = species_niche_pca[[sp]]$center[var],
                              scale = species_niche_pca[[sp]]$scale[var])
  }
  
  # Project these
  proj_cline = as.matrix(trans_cline) %*% species_niche_pca[[sp]]$rotation
  
  # Combine the projected and gbif data and do the hypervolume calcs
  combined_pc_res = rbind(proj_cline,
                          species_niche_pca[[sp]]$x) |>
    as.data.frame() |>
    mutate(group = c(rep('Sampled',nrow(proj_cline)),
                     rep('GBIF',nrow(species_niche_pca[[sp]]$x))))
  
  # # Hypervolume comparison by dynamic range boxes
  # overlap_res = dynRB_VPa(combined_pc_res[,c('group',colnames(combined_pc_res)[colnames(combined_pc_res) != 'group'])],
  #                         correlogram = F)
  # Hypervolume overlap of each of the PCs individually
  pc_overlap_res = dynRB_Pn(combined_pc_res[,c('group',colnames(combined_pc_res)[colnames(combined_pc_res) != 'group'])],)
  # Take a weighted mean, where weighting reflects the eigenvalues of the original PC
  pc_weights = species_niche_pca[[sp]]$sdev / sum(species_niche_pca[[sp]]$sdev)
  overlap_score = weighted.mean(as.numeric(pc_overlap_res$result[2,3:ncol(pc_overlap_res$result)]),
                                w = pc_weights)
  # # Return all the hypervolume overlap stats
  # out = data.table(species = sp,
  #                  prod_overlap = overlap_res$result$port_prod[2],
  #                  mean_overlap = overlap_res$result$port_mean[2])
  out = data.table(species = sp,
                   wmean_overlap = overlap_score)
  
  out
},mc.cores = 6) |>
  rbindlist()

# Save these resuts
write.csv(species_hypervolume_overlap,
          paste0('outputs/',output_name,'_gbif_niche_overlap.csv'),
          quote = F,row.names = F)
