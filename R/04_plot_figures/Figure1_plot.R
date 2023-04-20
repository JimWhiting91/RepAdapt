# Script for plotting all elements of figure 1
lib = c("ggridges","data.table","ggplot2","ggtree","ape","stars","sf","cowplot","pbmcapply")
sapply(lib,library,character.only=T)
source("R/repadapt_functions.R")

# General info
run_name = 230321
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"
OG_dir = "outputs/orthology/Results_221213_18_genomes_Ptaeda_isoforms_removed/"

# Fetch the OG_map
OG_map = data.table(readRDS(paste0("data/OF_OG_gea_gene_map_",output_name,".rds")))

# Plot a tree of all species
OG_tree = read.tree(paste0(OG_dir,"/Species_Tree/SpeciesTree_rooted.txt"))

# Image URLs for species...

# Fetch phylopic images...
genome_species_labs = c("Pseudotsuga menziesii",
                        "Picea abies\n(Picea obovata)\n(Picea glauca-engelmannii)",
                        "Pinus taeda\n(Pinus sylvestris)\n(Pinus contorta)",
                        "Panicum hallii",
                        "Amaranthus tuberculatus",
                        "Helianthus annuus\n(Helianthus argophyllus)\n(Helianthus petiolaris)",
                        "Medicago truncatula",
                        "Eucalyptus grandii*\n(Eucalyptus albens)\n(Eucalyptus sideroxylon)\n(Eucalyptus magnificata)",
                        "Quercus robur*\n(Quercus petraea)",
                        "Populus deltoides",
                        "Populus trichocarpa",
                        "Populus tremula",
                        "Arabis alpina",
                        "Capsella rubella",
                        "Boechera stricta",
                        "Arabidopsis thaliana",
                        "Arabidopsis halleri\n(Cardamine resedifolia)",
                        "Arabidopsis lyrata")
OG_tree$tip.label = genome_species_labs

# Remove conifer
OG_tree = drop.tip(OG_tree,tip = "Pseudotsuga menziesii")
# genome_species_url = list(Pseudotsuga_menziesii = c("https://upload.wikimedia.org/wikipedia/commons/thumb/d/d4/Pseudotsuga_menziesii_7456.JPG/800px-Pseudotsuga_menziesii_7456.JPG?20090212132550",
#                                                     "Walter Siegmund, CC BY-SA 3.0 <https://creativecommons.org/licenses/by-sa/3.0>, via Wikimedia Commons"),
#                           Picea_abies = c("https://upload.wikimedia.org/wikipedia/commons/9/9d/Kuusk_Keila-Paldiski_rdt_%C3%A4%C3%A4res.jpg",
#                                           "Ivar Leidus, CC BY-SA 3.0 EE <https://creativecommons.org/licenses/by-sa/3.0/ee/deed.en>, via Wikimedia Commons"),
#                           Pinus_taeda = c("https://upload.wikimedia.org/wikipedia/commons/e/ef/Pinus_taeda_crossett_exp_forest.jpg?20101215133757",
#                                           "U.S. Forest Service, Public domain, via Wikimedia Commons"),
#                           Panicum_hallii = c("http://inaturalist-open-data.s3.amazonaws.com/photos/14207623/medium.jpg",
#                                              '"Panicum hallii hallii" by samlutfy is licensed under CC BY 4.0. To view a copy of this license, visit https://creativecommons.org/licenses/by/4.0/?ref=openverse.'),
#                           Amaranthus_tuberculatus = c("https://upload.wikimedia.org/wikipedia/commons/e/e1/Amaranthus_tuberculatus_plant_%2802%29.jpg?20220909134643",
#                                                       "Jackson Campbell, CC BY-SA 4.0 <https://creativecommons.org/licenses/by-sa/4.0>, via Wikimedia Commons"),
#                           Helianthus_annuus = c("https://upload.wikimedia.org/wikipedia/commons/b/bd/Helianthus_annuus_exposed_2004-05-22.jpg?20051123172211",
#                                                 "Jon Sullivan, Public domain, via Wikimedia Commons"),
#                           Medicago_truncatula = c("https://upload.wikimedia.org/wikipedia/commons/thumb/5/5a/Medicago_truncatula_A17_branch.JPG/794px-Medicago_truncatula_A17_branch.JPG?20090723065509",
#                                                   "Ninjatacoshell, CC BY-SA 3.0 <https://creativecommons.org/licenses/by-sa/3.0>, via Wikimedia Commons"),
#                           Eucalyptus_grandis = c("https://upload.wikimedia.org/wikipedia/commons/d/d5/Eucalyptus_grandis_Kerewong_State_Forest_55_metres_tall.jpg?20110518062340",
#                                                  "Poyt448 Peter Woodard, CC BY-SA 3.0 <https://creativecommons.org/licenses/by-sa/3.0>, via Wikimedia Commons"),
#                           Quercus_petraea = c("https://upload.wikimedia.org/wikipedia/commons/thumb/8/89/Quercus_petraea_06.jpg/1216px-Quercus_petraea_06.jpg?20170317231252",
#                                               "Willow, CC BY-SA 3.0 <https://creativecommons.org/licenses/by-sa/3.0>, via Wikimedia Commons"),
#                           Populus_deltoides = c("https://upload.wikimedia.org/wikipedia/commons/thumb/1/11/Populus_deltoides_in_Golden_Valley_Tree_Park%2C_May_2022_02.jpg/900px-Populus_deltoides_in_Golden_Valley_Tree_Park%2C_May_2022_02.jpg?20220524060615",
#                                                 "Calistemon, CC BY-SA 4.0 <https://creativecommons.org/licenses/by-sa/4.0>, via Wikimedia Commons"),
#                           Populus_trichocarpa = c("https://upload.wikimedia.org/wikipedia/commons/thumb/3/39/Populus_trichocarpa_Umatilla.jpg/774px-Populus_trichocarpa_Umatilla.jpg?20090627194747",
#                                                   "Dave Powell, USDA Forest Service, CC BY 3.0 US <https://creativecommons.org/licenses/by/3.0/us/deed.en>, via Wikimedia Commons"),
#                           Populus_tremula = c("https://upload.wikimedia.org/wikipedia/commons/thumb/a/ac/PopulusTremula001.JPG/1600px-PopulusTremula001.JPG?20070620073531",
#                                               "I, Hugo.arg, CC BY-SA 3.0 <http://creativecommons.org/licenses/by-sa/3.0/>, via Wikimedia Commons"),
#                           Arabis_alpina = c("https://upload.wikimedia.org/wikipedia/commons/thumb/1/15/Arabis_alpina_kz10.jpg/1599px-Arabis_alpina_kz10.jpg?20200726133157",
#                                             "Krzysztof Ziarnek, Kenraiz, CC BY-SA 4.0 <https://creativecommons.org/licenses/by-sa/4.0>, via Wikimedia Commons"),
#                           Capsella_rubella = c("https://upload.wikimedia.org/wikipedia/commons/thumb/d/d2/Capsella_rubella_flower_%2818%29.jpg/793px-Capsella_rubella_flower_%2818%29.jpg?20211128181552",
#                                                "Jean-Jacques Houdré, CC BY-SA 2.0 FR <https://creativecommons.org/licenses/by-sa/2.0/fr/deed.en>, via Wikimedia Commons"),
#                           Boechera_stricta = c("https://live.staticflickr.com/3426/5840665045_9fea4b3c14_b.jpg",
#                                                '"Boechera stricta" by aspidoscelis is marked with CC0 1.0. To view the terms, visit https://creativecommons.org/publicdomain/zero/1.0/?ref=openverse.'),
#                           Arabidopsis_thaliana = c("https://upload.wikimedia.org/wikipedia/commons/thumb/0/03/Arabidopsis_thaliana_kz08.jpg/855px-Arabidopsis_thaliana_kz08.jpg?20220426171121",
#                                                    "Krzysztof Ziarnek, Kenraiz, CC BY-SA 4.0 <https://creativecommons.org/licenses/by-sa/4.0>, via Wikimedia Commons"),
#                           Arabidopsis_halleri = c("https://upload.wikimedia.org/wikipedia/commons/thumb/e/ea/Arabidopsis_halleri_flower_%2802%29.jpg/1210px-Arabidopsis_halleri_flower_%2802%29.jpg?20211104181227",
#                                                   "Andrea Moro, CC BY-SA 4.0 <https://creativecommons.org/licenses/by-sa/4.0>, via Wikimedia Commons"),
#                           Arabidopsis_lyrata = c("https://upload.wikimedia.org/wikipedia/commons/thumb/1/11/Arabidopsis_lyrata_-_Lyre_Leaf_Rock_Cress.jpg/450px-Arabidopsis_lyrata_-_Lyre_Leaf_Rock_Cress.jpg?20131214234811",
#                                                  "Fritzflohrreynolds, CC BY-SA 3.0 <https://creativecommons.org/licenses/by-sa/3.0>, via Wikimedia Commons")
# )
# phylopic_images <- ggimage::phylopic_uid(genome_species)
# phylopic_images2 <- phylopic_images
# phylopic_images$name = OG_tree$tip.label
# 
# p <- ggtree(OG_tree, branch.length = "none") %<+% phylopic_images
# 
# p + geom_tiplab(aes(label=label), offset = .2) + xlim(NA, 7) +
#   geom_tiplab(aes(image=uid), geom="phylopic",offset = 2.5)

# # imageInfo
# imageInfo = data.table(Newick_label = OG_tree$tip.label,
#                        vernacular_name = gsub("_"," ",names(genome_species_url)),
#                        imageURL = sapply(genome_species_url,'[[',1))

# p <- ggtree(OG_tree) %<+% imageInfo + xlim(NA, 6)
# p + geom_tiplab(aes(image = imageURL), geom="image", offset=rep(c(1,2),nrow(imageInfo)/2), align=T, size=.08, hjust=0) 
simple_OG_tree = ggtree(OG_tree) + geom_tiplab(size=12,geom="text",as_ylab = T)

# p + geom_tiplab(aes(label = vernacular_name),geom="label", offset=1, hjust=.5)


# library(ggimage)
# library(ggtree)
# url <- paste0("https://raw.githubusercontent.com/TreeViz/",
#               "metastyle/master/design/viz_targets_exercise/")
# x <- read.tree(paste0(url, "tree_boots.nwk"))
# info <- read.csv(paste0(url, "tip_data.csv"))

# Map elements ------------------------------------------------------------
# Fetch all of the sampling locations
climate_clines = list.files(paste0("outputs/GEA_res/"),recursive = T,pattern = "climate_cline")
climate_clines = grep("old|CoAdap",climate_clines,invert = T,value = T)
climate_clines = grep(run_name,climate_clines,value = T)

# Read in and fetch all sampling points
all_samplng_points = lapply(1:length(climate_clines),function(x){
  tmp = read.table(paste0("outputs/GEA_res/",climate_clines)[x],header = T)
  tmp$file = dirname(climate_clines[x])
  tmp
})
all_samplng_points_latlong = rbindlist(all_samplng_points)[,.(Lat,Long)]

# Plot world map for annual mean temp...
# Fetch relevant climate info
full_stars <- read_stars("data/worldclim/wc10/bio1.bil")
full_stars <- full_stars/10

# Plot global map...
temp_colors <- colorRampPalette(c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
global_temp_plot <- ggplot() + 
  geom_stars(data = full_stars) +
  scale_fill_gradientn(name = "Annual Mean Temp (BIOCLIM-1)",
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
        axis.ticks = element_blank())+
  geom_point(data=all_samplng_points_latlong,aes(x=Long,y=Lat),colour="black",size=1,alpha=0.8)

# Plot Japanese example of all of the clims...
japan_cline = read.table(paste0("outputs/GEA_res/",grep("ubota",climate_clines,value = T)),header = T)
set.seed(1000)
bioclim_plots = lapply(sort(sample(1:19,6)),function(x){
  print(x)
  
  # Fetch relevant climate info
  tmp_stars <- read_stars(paste0("data/worldclim/wc2-5/bio",x,".bil"))
  if(x < 12){
    tmp_stars <- tmp_stars/10
  }
  
  # Set the limits
  long_range <- max(japan_cline$Long) - min(japan_cline$Long)
  long_min <- min(japan_cline$Long) - 0.2*long_range
  long_max <- max(japan_cline$Long) + 0.2*long_range
  
  lat_range <- max(japan_cline$Lat) - min(japan_cline$Lat)
  lat_min <- min(japan_cline$Lat) - 0.2*lat_range
  lat_max <- max(japan_cline$Lat) + 0.2*lat_range
  
  # Filter the stars
  tmp_stars_part <- stars:::filter.stars(tmp_stars, x > long_min, x < long_max,y > lat_min,y < lat_max)
  
  # Plot
  climate_plot_tmp <- ggplot() + 
    geom_stars(data = tmp_stars_part) +
    scale_fill_gradientn(colors = temp_colors(5),
                         #limits = c(-7, 32),
                         na.value = "white") +
    coord_equal() +
    # scale_x_discrete(expand = c(0, 0)) +
    # scale_y_discrete(expand = c(0, 0)) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          axis.ticks = element_blank())+
    scale_x_continuous(limits = c(long_min,long_max))+
    scale_y_continuous(limits = c(lat_min,lat_max))+
    ggtitle(paste0("BIOCLIM-",x)) + 
    geom_point(data=japan_cline,aes(x=Long,y=Lat),colour="black",size=2,alpha=1)
  
  return(climate_plot_tmp)
})

bioclim_combined = cowplot::plot_grid(plotlist = bioclim_plots,
                                      nrow = 3)
# Combine bioclim and globe
cowplot::plot_grid(global_temp_plot,
                   bioclim_combined,
                   nrow = 2,rel_heights = c(1.2,1))

# prec clim change -----------------------------------------------
# Precipitation
raster_dir = "data/worldclim/wc2.1_2.5m_prec_1960-1969/"
prec_rasters1 = paste0(raster_dir,"/",list.files(raster_dir,pattern = "tif"))
raster_dir = "data/worldclim/wc2.1_2.5m_prec_2010-2018/"
prec_rasters2 = paste0(raster_dir,"/",list.files(raster_dir,pattern = "tif"))
# prec_rasters = c(prec_rasters,prec_rasters2)

sampling_data_list = list.files("data/VCFs/",pattern = "sampling_data.csv", recursive = TRUE,full.names = T)
sampling_data_kubota = grep('ubota',sampling_data_list,value = T)

long_lat_test = na.omit(unique(data.table(fread(sampling_data_kubota))[,.(Long,Lat)]))

# # Plot out some 1960s examples...
# plot_1960s = lapply(prec_rasters1[1:5],plot_climate_rasters,long_lat = long_lat_test)
# plot_2010s = lapply(prec_rasters2[1:5],plot_climate_rasters,long_lat = long_lat_test)
# cowplot::plot_grid(plotlist = plot_1960s)
# cowplot::plot_grid(plotlist = plot_2010s)

# Plot out some 1960s examples...
plot_1960s = plot_stacked_climate_rasters(prec_rasters1[1:12],long_lat = long_lat_test,stacking_val = 0.75) + 
  ggtitle(paste0("1961-1969 Precipitation\nN = ",length(prec_rasters1))) +
  theme(title = element_text(size = 12))
plot_2010s = plot_stacked_climate_rasters(prec_rasters2[1:12],long_lat = long_lat_test,stacking_val = 0.75) + 
  ggtitle(paste0("2010-2018 Precipitation\nN = ",length(prec_rasters2))) +
  theme(title = element_text(size = 12))

# # Plot these together as a test
# cowplot::plot_grid(plot_1960s,plot_2010s,
#                    ncol = 2)

#### Precipitation
prec_rasters = c(prec_rasters1,prec_rasters2)
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
prec_combine[grepl("201",prec_combine$date), date_group := "2010s"]

# For each population, calculate the effect size of change
prec_effectsize = rbindlist(lapply(unique(na.omit(prec_combine)[,pop]),function(x){
  tmp = prec_combine[pop == x,] %>% rstatix::wilcox_effsize(clim ~ date_group)
  data.table(pop = x,
             effsize = tmp$effsize)
}))
prec_effectsize = prec_effectsize[order(-effsize),]
prec_combine$pop_F = factor(prec_combine$pop,levels = rev(prec_effectsize$pop))
prec_effectsize$pop_F = factor(prec_effectsize$pop,levels = rev(prec_effectsize$pop))
# clim_change_example = ggplot(prec_combine,aes(pop_F,clim,fill = date_group)) +
#   geom_boxplot(outlier.colour = NA) +
#   scale_fill_brewer(palette = "Set1") +
#   theme_minimal() +
#   labs(y = "Annual Precipitation (mm)", x = "Sampling Location",fill = "") +
#   theme(panel.grid.major.x = element_blank(),
#         panel.grid.minor.x = element_blank(),
#         axis.text.x = element_blank(),
#         legend.position = c(0.95,0.95),
#         legend.direction = "horizontal",
#         legend.justification = c("right", "top"),
#         legend.box.background = element_rect(color="black", size=1),
#         legend.title = element_blank(),
#         legend.text = element_text(size = 12),
#         axis.text = element_text(size = 11),
#         axis.title = element_text(size = 14),
#         title = element_text(size = 14)) +
#   ggtitle("Example Precipitation Change 1960s-2010s")

clim_change_example = ggplot(prec_combine,aes(y = pop_F,x = clim,fill = date_group)) +
  stat_density_ridges(quantiles = 2,alpha = 0.5,quantile_lines = T) + 
  scale_fill_brewer(palette = "Set1") +
  theme_minimal() +
  labs(y = "Sampling Location", x = "Annual Precipitation (mm)",fill = "") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_blank(),
        # legend.position = c(0.95,0.95),
        legend.direction = "horizontal",
        # legend.justification = c("right", "top"),
        legend.box.background = element_rect(color="black", size=1,fill= 'white'),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 14),
        title = element_text(size = 12)) +
  ggtitle("Precipitation Change\n1960s-2010s")
clim_change_example

# Add effect size
clim_change_effect_example =ggplot(prec_effectsize,aes(y = pop_F,x = effsize)) +
  geom_bar(stat='identity',fill='gold2',colour = 'black') +
  theme_minimal() +
  labs(x = "Climate Change\nEffect Size",fill = "") +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 11),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_blank()) 

# Combine all of the climate change demo figs...
library(patchwork)
stacked_raster_patch = plot_1960s | plot_2010s
climchange_demo = plot_grid(stacked_raster_patch,
                            clim_change_example,
                            clim_change_effect_example,
                            align = 'h',axis = 'tblr',rel_widths = c(2,1,1),ncol = 3)

# clim_change_files = list.files("data/VCFs/",pattern = "clim", recursive = T)
# clim_change = rbindlist(lapply(paste0("data/VCFs/",clim_change_files),function(x){
#   tmp = fread(x)
#   tmp$dataset = dirname(x)
#   tmp
# }))

# Orthogroup Summary Figs -------------------------------------------------

# Remove the conifer...
OG_map = OG_map[genome != "Pmenziesii",]

# Count up numbers of paralogs in OG per genome
paralogN_stats = OG_map[,.(paralogN = nrow(.SD)),by = .(genome,Orthogroup)]
paralogN_stats$paralogN = as.character(paralogN_stats$paralogN)
paralogN_stats[!paralogN %in% c(as.character(1:10)),paralogN := ">10"]
total_paralogN_stats = paralogN_stats[,.(N = nrow(.SD)),by = .(genome,paralogN)]

# Produce a plot
OG_tree2 = read.tree(paste0(OG_dir,"/Species_Tree/SpeciesTree_rooted.txt"))
total_paralogN_stats$genome_F = factor(total_paralogN_stats$genome,levels = OG_tree2$tip.label)
total_paralogN_stats$paralogN_F = factor(total_paralogN_stats$paralogN,levels = c(">10",as.character(10:1)))
OG_paralogN_bars = ggplot(total_paralogN_stats,aes(y = genome_F,x = N,fill = paralogN_F)) +
  geom_bar(stat = "identity")+
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(x = "Number of Orthogroups (k)",fill = "Paralog N") + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom",
        # legend.direction = "vertical",
        # legend.justification = c("right", "top"),
        # legend.box.background = element_rect(color="black", size=1),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 14))+
  scale_x_continuous(labels = function(l) {trans = l / 1000;paste0(trans)}) +
  ggtitle("N Paralogs per Orthogroup")

# Now count up the number of genomes present in each orthogroup
OG_genome_counts = OG_map[,.(genomeN = length(unique(.SD$genome))),by = Orthogroup]
OG_genome_counts_merge = merge(OG_map,OG_genome_counts,by = "Orthogroup")
OG_genome_counts_merge = unique(OG_genome_counts_merge[,.(genome,Orthogroup,genomeN)])
total_OG_genome_counts = OG_genome_counts_merge[,.(N = nrow(.SD)),by = .(genome,genomeN)]
total_OG_genome_counts$genome_F = factor(total_OG_genome_counts$genome,levels = OG_tree2$tip.label)
total_OG_genome_counts$genomeN_F = factor(total_OG_genome_counts$genomeN,levels = as.character(1:17))
OG_genomeN_bars = ggplot(total_OG_genome_counts,aes(y = genome_F,x = N,fill = genomeN_F)) +
  geom_bar(stat = "identity")+
  scale_fill_viridis_d() +
  theme_minimal() +
  labs(x = "Number of Orthogroups (k)",fill = "Genome N") + 
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = "bottom",
        # legend.direction = "vertical",
        # legend.justification = c("right", "top"),
        # legend.box.background = element_rect(color="black", size=1),
        legend.text = element_text(size = 12),
        axis.text.x = element_text(size = 14),
        axis.title.x = element_text(size = 15),
        title = element_text(size = 14))+
  scale_x_continuous(labels = function(l) {trans = l / 1000;paste0(trans)}) +
  ggtitle("N Genomes per Orthogroup")


# Put everything together -------------------------------------------------
# pdf("figs/Figure1_study_OG_stats.pdf",width = 24,height = 11)
# cowplot::plot_grid(
#   cowplot::plot_grid(cowplot::plot_grid(global_temp_plot,bioclim_combined,ncol = 2,rel_widths = c(3,1),labels = c("A","B"),label_size = 32,label_x = c(0,-0.2)),
#                      # bioclim_combined,
#                      clim_change_example,
#                      nrow = 2,rel_heights = c(1.5,1),
#                      labels = c("","C"),label_size = 32,label_y = c(0,1.1)),
#   simple_OG_tree,
#   OG_paralogN_bars,
#   OG_genomeN_bars,
#   ncol = 4,align = "h",axis = 'tblr',labels = c("","D","E","F"),label_size = 32,label_x = c(0,0,-0.1,-0.1),
#   rel_widths = c(2.4,0.8,0.8,0.8)
# )
# dev.off()

pdf("figs/Figure1_study_OG_stats.pdf",width = 24,height = 11)
cowplot::plot_grid(
  cowplot::plot_grid(cowplot::plot_grid(global_temp_plot,bioclim_combined,ncol = 2,rel_widths = c(3,1),labels = c("A","B"),label_size = 32,label_x = c(0,-0.2)),
                     # bioclim_combined,
                     climchange_demo,
                     nrow = 2,rel_heights = c(1.5,1),
                     labels = c("","C"),label_size = 32,label_y = c(0,1.1)),
  simple_OG_tree,
  OG_paralogN_bars,
  OG_genomeN_bars,
  ncol = 4,align = "h",axis = 'tblr',labels = c("","D","E","F"),label_size = 32,label_x = c(0,0,-0.1,-0.1),
  rel_widths = c(2.4,0.8,0.8,0.8)
)
dev.off()



