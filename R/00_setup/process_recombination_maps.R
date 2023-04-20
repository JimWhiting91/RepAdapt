# Script for manipulating recombination maps and deriving rates
lib <- c("cowplot","poolr","mvmeta","qvalue","tidyr","ape","VGAM","ggExtra","pbmcapply","parallel","wCorr","data.table","ggplot2","viridis","ggridges","dplyr","readr")
sapply(lib,library,character.only=T)

# Functions
calc_rec_map_from_bad_map <- function(map_input){
  
  # Loop over chrs
  chr_out = rbindlist(lapply(unique(map_input[,1]),function(chr){
    # print(chr)
    
    # Subset
    map_tmp = map_input[map_input[,1] == chr,]
    
    # Mark out of order markers
    map_spline = smooth.spline(map_tmp[,c(2,3)],spar = 1)
    map_spline_dd = data.table(cm_pos = map_spline$y,
                               phys_pos = map_spline$x)
    
    
    # Find negative regions and just remove them
    rec_rates = rbindlist(lapply(2:nrow(map_spline_dd),function(x){
      
      out = data.table(chr = chr,
                       start = map_spline_dd[x-1,phys_pos],
                       end = map_spline_dd[x,phys_pos],
                       cM_start = map_spline_dd[x-1,cm_pos],
                       cM_end =  map_spline_dd[x,cm_pos])
      out$cM_Mb = (max(c(out$cM_end,out$cM_start)) - min(c(out$cM_end,out$cM_start)))/((out$end - out$start)/1000000)
      out
    }))
    
    rec_rates
  }))
  return(chr_out) 
}


# Aalpina -----------------------------------------------------------------
aalpina_map = read.csv("data/recombination_maps/A.alpina/A.alpina_genetic_map_v5.1.csv")

# Get rough rec rates...
aalpina_rec_rates = calc_rec_map_from_bad_map(aalpina_map[grep("chr",aalpina_map$Chromosome),c(1,2,5)])

# Check
ggplot(aalpina_rec_rates[cM_Mb < 50,],aes(start,cM_Mb))+geom_line()+
  facet_wrap(~chr,scales = "free")

# All look ok except for chr6 which is a mess...


# Alyrata -----------------------------------------------------------------
alyrata_map = read.csv("data/recombination_maps/A.lyrata/A.lyrara_hamala_2017_linkagemap.csv")

# We need to map the chroms
data.table(alyrata_map)[,.(len = max(mb)),by = lg][order(-len),]
alyrata_genome_chr_names = data.table(gea = c("Alyr_NW_003302555.1",
                                              "Alyr_NW_003302554.1",
                                              "Alyr_NW_003302553.1",
                                              "Alyr_NW_003302552.1",
                                              "Alyr_NW_003302551.1",
                                              "Alyr_NW_003302550.1",
                                              "Alyr_NW_003302549.1",
                                              "Alyr_NW_003302548.1"),
                                      rec = c(1:8))

# Get rough rec rates...
alyrata_rec_rates = calc_rec_map_from_bad_map(alyrata_map[,c(1,5,2)])

# Divide rec rates by 1 mil to account for double mb
alyrata_rec_rates$cM_Mb = alyrata_rec_rates$cM_Mb/1000000

# Check
ggplot(alyrata_rec_rates,aes(start,cM_Mb))+geom_line()+
  facet_wrap(~chr,scales = "free")

# Correct chrom names...
alyrata_rec_rates$chr = alyrata_genome_chr_names[match(alyrata_rec_rates$chr,alyrata_genome_chr_names$rec),gea]

# Fix cM to bp pos...
alyrata_rec_rates$start = alyrata_rec_rates$start * 1000000
alyrata_rec_rates$end = alyrata_rec_rates$end * 1000000


# A.thaliana --------------------------------------------------------------
athaliana_map = read.csv("data/recombination_maps/A.thaliana/A.thaliana_rowan_CO_map.csv")

# Get original chr names
athaliana_genome_chr_names = data.table(gea = c("Athal_NC_003070.9",
                                                "Athal_NC_003071.7",
                                                "Athal_NC_003074.8",
                                                "Athal_NC_003075.7",
                                                "Athal_NC_003076.8"),
                                        rec = 1:5)

# Split into 100kb windows and just count the number of crossovers...
athaliana_rec_rates = rbindlist(lapply(unique(athaliana_map$chr),function(chr){
  tmp = athaliana_map[athaliana_map$chr == chr,]
  
  # Split into 100kb intervals...
  tmp$wind = cut_interval(tmp$breakpoint.pos,length = 100000)
  
  # Crossovers per wind
  co_winds = data.frame(table(tmp$wind))
  
  # Tidy window
  co_winds$wind_start = as.integer(gsub("[","",gsub("(","",sapply(strsplit(as.character(co_winds$Var1),","),'[[',1),fixed = TRUE),fixed=TRUE))
  co_winds$wind_end = as.integer(gsub("]","",sapply(strsplit(as.character(co_winds$Var1),","),'[[',2),fixed = TRUE))
  
  # Save these...
  out = data.table(chr = tmp$chr[1],
                   start = co_winds$wind_start,
                   end = co_winds$wind_end,
                   co_100k = co_winds$Freq)
}))

# Check
ggplot(athaliana_rec_rates,aes(start,co_100k))+
  geom_line()+
  facet_wrap(~chr,scales = "free_x")

# Make (wrongly) named column
athaliana_rec_rates$cM_Mb <- athaliana_rec_rates$co_100k

# Correct chrom names...
athaliana_rec_rates$chr = athaliana_genome_chr_names[match(athaliana_rec_rates$chr,athaliana_genome_chr_names$rec),gea]


# H.annuus ----------------------------------------------------------------
hannuus_map = fread("data/recombination_maps/H.annuus/Ha412HOv2.0-20181130.Nov22k22.geneticmap.extradivisions.txt")

ggplot(hannuus_map,aes(pos,cM))+
  geom_line()+
  facet_wrap(~chr,scales = "free")

# Get rough rec rates...
hannuus_rec_rates = calc_rec_map_from_bad_map(data.frame(hannuus_map))

# Check
ggplot(hannuus_rec_rates,aes(start,cM_Mb))+geom_line()+
  facet_wrap(~chr,scales = "free")

# # Find negative regions and just remove them
# hannuus_rec_rates = rbindlist(lapply(unique(hannuus_map$chr),function(chr){
#   map_tmp = hannuus_map[chr == chr,]
#   
#   rec_rates = rbindlist(lapply(2:nrow(map_tmp),function(x){
#     
#     out = data.table(chr = chr,
#                      start = map_tmp[x-1,pos],
#                      end = map_tmp[x,pos],
#                      cM_start = map_tmp[x-1,cM],
#                      cM_end =  map_tmp[x,cM])
#     out$cM_Mb = (max(c(out$cM_end,out$cM_start)) - min(c(out$cM_end,out$cM_start)))/((out$end - out$start)/1000000)
#     out
#   }))
#   
# }))

# Save all the recombination rates as an R object
out_list = list(Aalpina = aalpina_rec_rates,
                Alyrata = alyrata_rec_rates,
                Athaliana = athaliana_rec_rates,
                Hannuus = hannuus_rec_rates)
saveRDS(out_list,"outputs/processed_recombination_rates.rds")

