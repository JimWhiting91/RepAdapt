# General metadata --------------------------------------------------------
species_genome_map = data.table(species = c("Pseudotsuga menziesii",
                                            "Picea abies","Picea obovata","Picea glaucaxengelmannii",
                                            "Pinus sylvestris","Pinus contorta",
                                            "Panicum hallii",
                                            "Amaranthus tuberculatus",
                                            "Helianthus annuus","Helianthus argophyllus","Helianthus petiolaris",
                                            "Medicago truncatula",
                                            "Eucalyptus albens","Eucalyptus sideroxylon","Eucalyptus magnificata",
                                            "Quercus petraea",
                                            "Populus deltoides",
                                            "Populus trichocarpa",
                                            "Populus tremula",
                                            "Arabis alpina",
                                            "Capsella rubella","Capsella bursa-pastoris",
                                            "Boechera stricta",
                                            "Arabidopsis thaliana",
                                            "Arabidopsis halleri","Cardamine resedifolia",
                                            "Arabidopsis lyrata"),
                                genome = c("Pmenziesii",
                                           "Pabies","Pabies","Pabies",
                                           "Ptaeda","Ptaeda",
                                           "Phallii",
                                           "Atuberculatus",
                                           "Hannuus","Hannuus","Hannuus",
                                           "Mtruncatula",
                                           "Egrandis","Egrandis","Egrandis",
                                           "Qpetraea",
                                           "Pdeltoides",
                                           "Ptrichocarpa",
                                           "Ptremula",
                                           "Aalpina",
                                           "Crubella","Crubella",
                                           "Bstricta",
                                           "Athaliana",
                                           "Ahalleri","Ahalleri",
                                           "Alyrata"))


#### RepAdapt Function Library

# Function takes some data and condenses down to Orthogroup-level empirical pvals
condenseToOGPvals <- function(data,OG_col,pval_col){
  
  data$Orthogroup = data.frame(data)[,OG_col]
  data$tmp_p = data.frame(data)[,pval_col]
  
  # Group summarise
  OG_minP = data[,.(minP = min(.SD$tmp_p),
                    N = nrow(.SD)),by = Orthogroup]
  
  # Convert to DS corrected
  OG_minP$minP_DS =  1 - (1 - OG_minP$minP)^OG_minP$N
  
  # Add epvalue
  OG_minP$minP_DS_epvalue = empPvals(-OG_minP$minP_DS,-OG_minP$minP_DS)
  OG_minP
}

# Function extracts climate values for a given lat long from a raster
extract_climate_values = function(long_lat,raster_path){
  library(raster)
  library(sf)
  clim = raster(raster_path)
  
  # Match projections and extract clim values
  coordinates(long_lat) <- ~Long+Lat
  mypoints = SpatialPoints(long_lat,proj4string = crs("+init=epsg:4326"))
  myproj = crs(clim)
  points.proj = spTransform(mypoints, myproj)
  raster::extract(clim, points.proj)
}

#' Rotate simple features for 3D layers
#' Rotates a simple features layer using a shear matrix transformation on the 
#' \code{geometry} column. This can get nice for visualisation and works with
#' points, lines and polygons.
#'
#' @param data an object of class \code{sf}
#' @param x_add integer; x value to move geometry in space
#' @param y_add integer; x value to move geometry in space
#'
#' #' @importFrom magrittr %>%

rotate_sf <- function(data, x_add = 0, y_add = 0) {
  
  shear_matrix <- function (x) { 
    matrix(c(2, 1.2, 0, 1), 2, 2) 
  }
  
  rotate_matrix <- function(x) { 
    matrix(c(cos(x), sin(x), -sin(x), cos(x)), 2, 2) 
  }
  
  data %>% 
    dplyr::mutate(
      geometry = 
        .$geometry * shear_matrix() * rotate_matrix(pi / 20) + c(x_add, y_add)
    )
}

# Plot rasters based on input and long_lat info
# Stack these if there are multiple rasters...
plot_climate_rasters = function(long_lat,raster_path){
  temp_colors <- colorRampPalette(c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
  
  library(raster)
  library(sf)
  library(stars)
  # clim = raster(raster_path)
  # extent(clim) = extent(long_lat)
  # clim_trim = setExtent(clim,extent(long_lat))
  
  long_min = min(long_lat$Long)
  long_max = max(long_lat$Long)
  lat_min = min(long_lat$Lat)
  lat_max = max(long_lat$Lat)
  long_extra = (long_max - long_min) * 0.1
  lat_extra = (lat_max - lat_min) * 0.1
  
  clim = stars:::read_stars(raster_path)
  tmp_stars_part <- stars:::filter.stars(clim, 
                                         x > long_min - long_extra, x < long_max + long_extra,
                                         y > lat_min - lat_extra, y < lat_max + lat_extra)
  
  ggplot() + 
    geom_stars(data = tmp_stars_part) +
    scale_fill_gradientn(colors = temp_colors(5),
                         #limits = c(-7, 32),
                         na.value = "white") +
    coord_equal() +
    # scale_x_discrete(expand = c(0, 0)) +
    # scale_y_discrete(expand = c(0, 0)) +
    theme_bw() +
    # theme_minimal() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          panel.border = element_rect(size = 2),
          plot.margin=grid::unit(c(0,0,0,0), "mm"))+
    scale_x_continuous(limits = c(long_min - long_extra,long_max + long_extra))+
    scale_y_continuous(limits = c(lat_min - lat_extra,lat_max + lat_extra))+
    geom_point(data=long_lat,aes(x=Long,y=Lat),colour="black",size=2,alpha=1) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0))
}

# This takes a set of raster paths and plots them as a stack...
plot_stacked_climate_rasters = function(long_lat,
                                        raster_path,
                                        stacking_val = 1){
  
  # Use these colours
  temp_colors <- colorRampPalette(c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
  
  library(raster)
  library(sf)
  library(stars)
  # clim = raster(raster_path)
  # extent(clim) = extent(long_lat)
  # clim_trim = setExtent(clim,extent(long_lat))
  
  long_min = min(long_lat$Long)
  long_max = max(long_lat$Long)
  lat_min = min(long_lat$Lat)
  lat_max = max(long_lat$Lat)
  long_extra = (long_max - long_min) * 0.1
  lat_extra = (lat_max - lat_min) * 0.1
  
  # Read in the first raster, this will go at the top
  clim = stars:::read_stars(raster_path)
  stars_part <- stars:::filter.stars(clim, 
                                     x > long_min - long_extra, x < long_max + long_extra,
                                     y > lat_min - lat_extra, y < lat_max + lat_extra)
  
  # Transform to sf and tilt those mfs
  sm <- matrix(c(2, 1.2, 0, 1), 2, 2)
  clim_tilt <- st_as_sf(stars_part)
  clim_tilt$geometry <- clim_tilt$geometry * sm
  st_crs(clim_tilt) <- st_crs(st_as_sf(stars_part))
  
  # Save a list of rasters but extract the top one and rename
  names(clim_tilt)[1] = 'climate'
  all_clim_tilt = clim_tilt
  clim_tilt = clim_tilt |> dplyr::select("climate")
  
  # Read in all the rasters...
  if(length(raster_path) > 1){
    other_rasters = lapply(2:(length(all_clim_tilt) - 1),function(i){
      # print(i)
      # Subset...
      clim_tilt_tmp = all_clim_tilt |> dplyr::select(all_of(names(all_clim_tilt)[i]))
      clim_tilt_tmp$geometry <- clim_tilt_tmp$geometry + c(0, -1 * stacking_val * (i - 1))
      st_crs(clim_tilt_tmp) <- st_crs(st_as_sf(clim_tilt))
      names(clim_tilt_tmp)[1] = 'climate'
      return(clim_tilt_tmp)
    })
  }
  
  # # Read in all the rasters...
  # if(length(raster_path) > 1){
  #   other_rasters = lapply(2:length(raster_path),function(i){
  #     
  #     # Read in
  #     clim_tmp = stars:::read_stars(raster_path[i])
  #     tmp_stars_part <- stars:::filter.stars(clim_tmp, 
  #                                            x > long_min - long_extra, x < long_max + long_extra,
  #                                            y > lat_min - lat_extra, y < lat_max + lat_extra)
  #     
  #     # Tilt the same way
  #     clim_tilt_tmp <- st_as_sf(tmp_stars_part)
  #     clim_tilt_tmp$geometry <- clim_tilt_tmp$geometry * sm + c(0, -1 * stacking_val * (i - 1))
  #     st_crs(clim_tilt_tmp) <- st_crs(st_as_sf(clim_tilt))
  #     names(clim_tilt_tmp)[1] = 'climate'
  #     return(clim_tilt_tmp)
  #   })
  # }
  
  # # pl
  # clim_tilt2 = clim_tilt
  # st_geometry(clim_tilt2) <- st_geometry(clim_tilt2) + c(0, -1)
  # st_crs(clim_tilt2) <- st_crs(clim_tilt)
  
  # Plot the first one...
  p1 = ggplot() + 
    geom_sf(data = clim_tilt,size = 10,colour = 'black',lwd = 1.5) +
    geom_sf(data = clim_tilt,aes(fill = climate,
                                 colour = climate)) +
    scale_fill_gradientn(colors = temp_colors(5),
                           #limits = c(-7, 32),
                           na.value = "white") + 
    scale_colour_gradientn(colors = temp_colors(5),
                         #limits = c(-7, 32),
                         na.value = "white")
  # Now add each of these in reverse order so they stack...
  for(i in length(other_rasters):1){
    tmp_tilt = other_rasters[[i]]
    p1 = p1 + geom_sf(data = tmp_tilt,lwd = 1.5,colour = 'black') +
      geom_sf(data = tmp_tilt,aes(colour = climate,
                                  fill = climate)) 
  }
  
  # Add top layer back on
  p1 = p1 + 
    geom_sf(data = clim_tilt,size = 10,colour = 'black',lwd = 1.5) +
    geom_sf(data = clim_tilt,aes(fill = climate,
                                 colour = climate))
  
  # Get our points ready...
  long_lat_sf = st_as_sf(long_lat, coords = c(1:2))
  long_lat_sf$geometry = long_lat_sf$geometry * sm
  st_crs(long_lat_sf) <- st_crs(st_as_sf(clim_tilt))
  

  # Add on original points and do all remaining formatting
  p1 + theme_minimal() +
    # theme_minimal() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    geom_sf(data = long_lat_sf,colour = 'black',size = 3) +
    geom_sf(data = long_lat_sf,colour = 'white',size = 2)

}


# Main function to do prec tmax calcs
calculate_climchange_from_longlat = function(sampling_path){
  
  # Fetch sampling data from kubota
  long_lat_test = na.omit(unique(data.table(fread(sampling_path))[,.(Long,Lat)]))
  
  # Fetch all the data through time...
  #### Precipitation
  prec_time_clim = pbmclapply(prec_rasters,function(raster_in){
    clim_tmp = extract_climate_values(raster_path = raster_in,
                                      long_lat = long_lat_test)
    year_tmp =  stringr::str_split(raster_in, "_|\\.")[[1]]
    out = data.table(date = year_tmp[length(year_tmp) - 1],
                     clim = clim_tmp,
                     pop = apply(long_lat_test,1,paste,collapse = "_"))
  },mc.cores = 4)
  
  prec_combine = rbindlist(prec_time_clim)
  prec_combine$date_F = factor(prec_combine$date,levels = unique(sort(prec_combine$date)))
  
  # Set before and after
  prec_combine$date_group = "1960s"
  prec_combine[grepl("201",prec_combine$date),date_group := "2010s"]
  
  # For each population, calculate the effect size of change
  prec_effectsize = rbindlist(lapply(unique(na.omit(prec_combine)[,pop]),function(x){
    # print(x)
    tmp = prec_combine[pop == x,] %>% rstatix::wilcox_effsize(clim ~ date_group)
    data.table(pop = x,
               effsize = tmp$effsize)
  }))
  
  #### Max Temp
  tmax_time_clim = pbmclapply(tmax_rasters,function(raster_in){
    clim_tmp = extract_climate_values(raster_path = raster_in,
                                      long_lat = long_lat_test)
    year_tmp =  stringr::str_split(raster_in, "_|\\.")[[1]]
    out = data.table(date = year_tmp[length(year_tmp) - 1],
                     clim = clim_tmp,
                     pop = apply(long_lat_test,1,paste,collapse = "_"))
  },mc.cores = 4)
  
  tmax_combine = rbindlist(tmax_time_clim)
  tmax_combine$date_F = factor(tmax_combine$date,levels = unique(sort(tmax_combine$date)))
  
  # Set before and after
  tmax_combine$date_group = "1960s"
  tmax_combine[grepl("201",tmax_combine$date), date_group := "2010s"]
  
  # For each population, calculate the effect size of change
  tmax_effectsize = rbindlist(lapply(unique(na.omit(tmax_combine)[,pop]),function(x){
    tmp = tmax_combine[pop == x,] %>% rstatix::wilcox_effsize(clim ~ date_group)
    data.table(pop = x,
               effsize = tmp$effsize)
  }))
  
  ### Merge all together
  clim_change_merge = merge(prec_effectsize,tmax_effectsize, by = "pop")
  colnames(clim_change_merge) = c("pop","prec_clim_change","tmax_clim_change")
  return(clim_change_merge)
}


# Get the species order from the species tree... --------------------------
get_species_order = function(){
  OF_tree = read.tree("outputs/orthology/Results_221213_18_genomes_Ptaeda_isoforms_removed/Species_Tree/SpeciesTree_rooted.txt")
  # Fetch our genome dataset map
  dataset_meta <- read.table("metadata/vcf_genome_gff_220830_map.txt")
  dataset_meta <- na.omit(dataset_meta[,c(8,5)])
  colnames(dataset_meta) <- c("dataset","genome")
  
  # Reorder metadata...
  dataset_meta_ordered <- NULL
  for(genome in OF_tree$tip.label){
    dataset_meta_ordered <- rbind(dataset_meta_ordered,dataset_meta[dataset_meta$genome == genome,])
  }
  
  # Reduce down for species...
  dataset_meta_ordered$species = paste0(sapply(strsplit(dataset_meta_ordered$dataset,'_'),'[[',1),
                                        " ",
                                        sapply(strsplit(dataset_meta_ordered$dataset,'_'),'[[',2))
  
  
  species_meta_ordered = unique(dataset_meta_ordered[,c("genome","species")])
  return(species_meta_ordered)
}

# moran_OG_tree <- function(OG,OG_tree_dir,climate_var,plot.tree=T){
#   
#   # Fetch gea_res
#   all_gea_res <- all_picmin_res[[climate_var]]$OG_maxP
#   
#   # Fetch the tree
#   if(file.exists(paste0(OG_tree_dir,"/",OG,"_tree.txt"))){
#     tmp_tree <- read.tree(paste0(OG_tree_dir,"/",OG,"_tree.txt"))
#     
#     # Build tree metadata
#     OG_tree_metadata <- data.frame(gene_name=tmp_tree$tip.label,
#                                    tip.pos=1:length(tmp_tree$tip.label),
#                                    OG=OG)
#     
#     # Add on OF codes
#     OG_tree_metadata <- merge(OG_tree_metadata,OG_map[,c("OF_gene","gene_name")],by="gene_name")
#     OG_tree_metadata <- OG_tree_metadata[order(OG_tree_metadata$tip.pos),]
#     
#     # Fetch all of the relevant GEA results...
#     OG_GEA_res <- rbindlist(lapply(unique(all_gea_res$dataset),function(gea){
#       
#       # Subset for OG genes
#       gea_sub <- gea[full_OF %in% OG_tree_metadata$OF_gene,]
#       
#       if(nrow(gea_sub)==0){
#         return(NULL)
#       } else {
#         gea_sub$OF_gene <- gea_sub$full_OF
#         gea_sub <- merge(gea_sub,OG_tree_metadata[,c("tip.pos","OF_gene")],by="OF_gene")
#         return(gea_sub)
#       }
#       
#       
#     }))
#     
#     # Get branch distances
#     OG_branch_dists <- cophenetic.phylo(tmp_tree)
#     rownames(OG_branch_dists) <- OG_tree_metadata$OF_gene[match(rownames(OG_branch_dists),OG_tree_metadata$gene_name)]
#     colnames(OG_branch_dists) <- OG_tree_metadata$OF_gene[match(colnames(OG_branch_dists),OG_tree_metadata$gene_name)]
#     
#     # Only keep tips that represent the minimum p value from each dataset
#     tree_min_tips <- OG_GEA_res[,.(minP=min(.SD$pvalue),
#                                    OF_gene=(.SD$OF_gene[which(.SD$pvalue == min(.SD$pvalue))]),
#                                    tip.pos=(.SD$tip.pos[which(.SD$pvalue == min(.SD$pvalue))])),by=dataset]
#     
#     if(plot.tree){
#       
#       # Reduce tree
#       
#       gea_mat <- matrix(nrow = length(unique(tree_min_tips$tip.pos)),ncol=max(table(tree_min_tips$tip.pos)))
#       tmp_tree2 = tmp_tree
#       tmp_tree2$tip.label <- OG_tree_metadata$OF_gene
#       tmp_tree_reduced <- drop.tip(tmp_tree2,tip = tmp_tree2$tip.label[!(tmp_tree2$tip.label %in% tree_min_tips$OF_gene)])
#       OG_tree_metadata_reduced <- data.frame(gene_name=tmp_tree_reduced$tip.label,
#                                              tip.pos=1:length(tmp_tree_reduced$tip.label),
#                                              OG=OG)
#       
#       # plot(tmp_tree_reduced)
#       
#       # Make heatmap
#       for(i in 1:nrow(OG_tree_metadata_reduced)){
#         to_fill = table(tree_min_tips$OF_gene)[OG_tree_metadata_reduced$gene_name[i]]
#         for(j in 1:to_fill){
#           gea_mat[i,j] <- tree_min_tips[tree_min_tips$OF_gene == OG_tree_metadata_reduced$gene_name[i],]$minP[j]
#         }
#       }
#       rownames(gea_mat) <- names(table(tree_min_tips$OF_gene)[OG_tree_metadata_reduced$gene_name])
#       # gea_melt <- suppressWarnings(reshape2::melt(gea_mat))
#       # gea_heatmap <- ggplot(na.omit(gea_melt),(aes(y=Var1,x=Var2,fill=-log10(value))))+
#       #   geom_tile(na.rm = T)+
#       #   scale_fill_viridis(option="A")+
#       #   theme_void()+
#       #   theme(legend.title = element_text(size=18),
#       #         legend.text = element_text(size=14))+
#       #   labs(fill="GEA\nPval")
#       
#       # Make tree with heatmap...
#       gene_tree = ggtree(tmp_tree_reduced, branch.length='none')+ 
#         theme_void()
#       out_fig <- gheatmap(gene_tree, as.data.frame(-log10(gea_mat)),
#                           colnames=F, legend_title="GEA Pval")+
#         scale_fill_viridis(option="A")
#       
#     }
#     
#     # Make the new weights matrix
#     weight_mat <- matrix(ncol=length(tree_min_tips$minP),nrow=length(tree_min_tips$minP))
#     for(i in 1:nrow(weight_mat)){
#       for(j in 1:ncol(weight_mat)){
#         weight_mat[i,j] <- OG_branch_dists[tree_min_tips$OF_gene[i],tree_min_tips$OF_gene[j]]
#       }
#     }
#     
#     # Now calculate moran's I...
#     test_vector <- -log10(tree_min_tips$minP)
#     names(test_vector) <- tree_min_tips$OF_gene
#     OG_branch_weights <- 1/OG_branch_dists[names(test_vector),names(test_vector)]
#     OG_branch_weights[!is.finite(OG_branch_weights)] <- 0
#     if(sum(OG_branch_weights) > 0){
#       out <- data.frame(Moran.I(test_vector,weight = OG_branch_weights,scaled = T))
#       out$climate_var = climate_var
#       out$OG = OG
#       out$n_tips = length(test_vector)
#       
#       
#       if(plot.tree){
#         return(list(res=out,tree_fig=out_fig))
#       } else {
#         return(out)
#       }
#     } else {
#       return(NULL) 
#     }
#   } else {
#     return(NULL)
#   }
# }

# Function Library --------------------------------------------------------
mirror_wza_sd_calc <- function(wza_vector,snp_count_vector){
  
  # First, we take the lower half of the wza vector and mirror it...
  wza_maxy <- which.max(density(wza_vector)$y) # Find the density peak
  wza_density_max <- density(wza_vector)$x[wza_maxy] # Find the wza at the density
  to_keep <- wza_vector < wza_density_max
  lowerhalf_wza_std <- wza_vector[to_keep] - max(wza_vector[to_keep])
  lowerhalf_wza_mirror <- c(lowerhalf_wza_std,-1*lowerhalf_wza_std) # Make a mirror half
  
  # Now return a data.table with mean snp count and
  out = data.table(wza_sd = sd(lowerhalf_wza_mirror),
                   snp_count = mean(snp_count_vector[to_keep]))
}

GO_enrichment_orthogroups = function(OG_input,GO_universe,n_cores = 1){
  
  # Fetch the go terms...
  to_test = unique(GO_universe[Orthogroup %in% OG_input,go_id])
  to_test = to_test[to_test != ""]
  
  # Count totals in the universe
  go_totals = data.table(table(GO_universe[go_id %in% to_test,go_id]))
  
  # Only test where we have at least 2 to test...
  go_totals_in_set = table(GO_universe[Orthogroup %in% OG_input,go_id])
  to_test = to_test[to_test %in% names(go_totals_in_set[go_totals_in_set > 1])]
  
  # Hypergeometric test for each GO...
  hyper_res = rbindlist(pbmclapply(to_test,function(tmp_go){
    nhits = nrow(GO_universe[go_id == tmp_go & Orthogroup %in% OG_input,])
    whiteball_n = length(unique(GO_universe[go_id == tmp_go,Orthogroup]))
    blackball_n = sum(!(unique(GO_universe$Orthogroup) %in% unique(GO_universe[go_id == tmp_go,Orthogroup])))
    ndraws = length(OG_input)
    
    # Test with fishers exact test...
    test_mat <- data.frame(
      "convergent" = c(nhits, ndraws - nhits),
      "total" = c(whiteball_n, blackball_n),
      row.names = c("with", "without"),
      stringsAsFactors = FALSE
    )
    
    data.table(go_id = tmp_go,
               nhits = nhits,
               total = whiteball_n,
               obs = round(nhits/ndraws,3),
               exp = round(whiteball_n/blackball_n,3),
               pval = fisher.test(test_mat)$p)
  },mc.cores = n_cores))
  
  # Add FDR
  hyper_res$fdr = p.adjust(hyper_res$pval,method = "fdr")
  
  # Merge back with descriptions
  hyper_res = merge(hyper_res,unique(GO_universe[,.(go_id,go_name)]))[order(fdr),]
  
  hyper_res[fdr < 0.1,]
  
}

# Repeatability Testing ---------------------------------------------------
DunnSidak_pvals = function(pval_vector){
  1 - (1 - min(pval_vector))^length(pval_vector)
}

stouffer_pvals = function(pval_vector){
  sum(qnorm(1 - pval_vector)) / sqrt(length(pval_vector))
}

stouffer_pvals_noMin = function(pval_vector){
  pval_vector2 = pval_vector[-which.min(pval_vector)]
  sum(qnorm(1 - pval_vector2)) / sqrt(length(pval_vector2))
}



