# Script for plotting all elements of figure 1
lib = c("ggridges","data.table","ggplot2","ggtree","ape","stars","sf","cowplot","pbmcapply")
sapply(lib,library,character.only=T)
source("R/repadapt_functions.R")

# General info
run_name = 230321
output_name = "25species_fixedAlyrataPabiesPobovata_OFcodes"
OG_dir = "outputs/orthology/Results_221213_18_genomes_Ptaeda_isoforms_removed//"

# Fetch the OG_map
OG_map = data.table(readRDS(paste0("data/OF_OG_gea_gene_map_",output_name,".rds")))

# Plot a tree of all species
OG_tree = read.tree(paste0(OG_dir,"/Species_Tree/SpeciesTree_rooted.txt"))

# Fetch phylopic images...
genome_species_labs = c(Pmenziesii = "Pseudotsuga menziesii",
                        Pabies = "Picea abies\n(Picea obovata)\n(Picea glauca-engelmannii)",
                        Ptaeda = "Pinus taeda\n(Pinus sylvestris)\n(Pinus contorta)",
                        Phallii = "Panicum hallii",
                        Atubercatus = "Amaranthus tuberculatus",
                        Hannuus = "Helianthus annuus\n(Helianthus argophyllus)\n(Helianthus petiolaris)",
                        Mtruncatula = "Medicago truncatula",
                        Egrandis = "Eucalyptus grandis*\n(Eucalyptus albens)\n(Eucalyptus sideroxylon)\n(Eucalyptus magnificata)",
                        Qpetraea = "Quercus robur*\n(Quercus petraea)",
                        Pdeltoides = "Populus deltoides",
                        Ptrichocarpa = "Populus trichocarpa",
                        Ptremula = "Populus tremula",
                        Aalpina = "Arabis alpina",
                        Crubella = "Capsella rubella",
                        Bstricta = "Boechera stricta",
                        Athaliana = "Arabidopsis thaliana",
                        Ahalleri = "Arabidopsis halleri\n(Cardamine resedifolia)",
                        Alyrata = "Arabidopsis lyrata")
for(i in 1:length(OG_tree$tip.label)){
  OG_tree$tip.label[i] = genome_species_labs[OG_tree$tip.label[i]]
}

# Remove conifer
OG_tree = drop.tip(OG_tree,tip = "Pseudotsuga menziesii")

# Add some node labels
OG_tree$node.label = ''
OG_tree$node.label[1] = '330.3\nmya'
OG_tree$node.label[2] = '130.2 mya'
OG_tree$node.label[3] = '159.6 mya'
OG_tree$node.label[4] = '118 mya'
OG_tree$node.label[15] = '15.97 mya'
OG_tree$node.label[7] = '25.97 mya'


simple_OG_tree = ggtree(OG_tree) + 
  geom_nodepoint() +
  geom_label2(aes(subset = node == 18,label = '330.3\nMYA'),fill = 'orange3',alpha = 0.5,nudge_x = 0.05) +
  geom_label2(aes(subset = node == 19,label = '130.2\nMYA'),fill = 'orange3',alpha = 0.5) +
  geom_label2(aes(subset = node == 20,label = '159.6\nMYA'),fill = 'orange3',alpha = 0.5) +
  geom_label2(aes(subset = node == 32,label = '15.97\nMYA'),fill = 'orange3',alpha = 0.5) +
  geom_label2(aes(subset = node == 24,label = '25.97\nMYA'),fill = 'orange3',alpha = 0.5) +
  geom_tiplab(size=12,geom="text",as_ylab = T)


simple_OG_tree

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

sampling_data_list = list.files("data/VCFs/",pattern = "sampling_data.csv", recursive = TRUE,full.names = T)
sampling_data_kubota = grep('ubota',sampling_data_list,value = T)

long_lat_test = na.omit(unique(data.table(fread(sampling_data_kubota))[,.(Long,Lat)]))


# Plot out some 1960s examples...
plot_1960s = plot_stacked_climate_rasters(prec_rasters1[1:12],long_lat = long_lat_test,stacking_val = 0.75) + 
  ggtitle(paste0("1961-1969 Precipitation\nN = ",length(prec_rasters1))) +
  theme(title = element_text(size = 12))
plot_2010s = plot_stacked_climate_rasters(prec_rasters2[1:12],long_lat = long_lat_test,stacking_val = 0.75) + 
  ggtitle(paste0("2010-2018 Precipitation\nN = ",length(prec_rasters2))) +
  theme(title = element_text(size = 12))


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

# Plot a demonstration of how climate change variables look
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
climchange_demo = cowplot::plot_grid(stacked_raster_patch,
                                     clim_change_example,
                                     clim_change_effect_example,
                                     align = 'h',axis = 'tblr',rel_widths = c(2,1,1),ncol = 3)

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

pdf("figs/Figure1_study_OG_stats.pdf",width = 18,height = 22)
cowplot::plot_grid(
  cowplot::plot_grid(cowplot::plot_grid(global_temp_plot,bioclim_combined,ncol = 2,rel_widths = c(3,1),labels = c("A","B"),label_size = 32,label_x = c(0,-0.2)),
                     # bioclim_combined,
                     climchange_demo,
                     nrow = 2,rel_heights = c(2,1),
                     labels = c("","C"),label_size = 32,label_y = c(0,1.1)),
  cowplot::plot_grid(simple_OG_tree,
                     OG_paralogN_bars,
                     OG_genomeN_bars,
                     ncol = 3,align = "h",axis = 'tblr',labels = c("D","E","F"),label_size = 32,label_x = c(0,-0.1,-0.1)),
  nrow = 2,rel_heights = c(1.2,1)
)
dev.off()



