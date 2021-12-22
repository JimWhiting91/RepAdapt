#### HMM classification of peaks into clusters...
lib <- c("depmixS4","R.utils","doParallel","pbapply","vegan","data.table","ggplot2","dplyr","tidyr")
dummy <- suppressPackageStartupMessages(lapply(lib,function(x){ 
  if (!require(x,character.only = T)){ 
    install.packages(x,repos='https://utstat.toronto.edu/cran/')
    library(x,character.only = T)
  } else {
    library(x,character.only = T)
  }
}))

# Fetch command args ------------------------------------------------------
args <- commandArgs(asValues = T,excludeReserved = T)[-1]
print(args)

res_dir <- paste0("outputs/GEA_res/",args$dataset_dir)
n_cores <- args$n_cores
window_sizes <- as.integer(unlist(strsplit(args$window_sizes,",")))
n_states <- as.integer(unlist(strsplit(args$n_states,",")))

if(any(is.null(c(res_dir)))){
  stop("ERROR: No GEA results provided")
}
if(is.null(n_cores)){
  n_cores <- 1
}
if(is.null(window_sizes)){
  window_sizes <- 100
}
if(is.null(n_states)){
  n_states <- c(3,5)
}

# Set up parallel cores
cl <- makeCluster(n_cores, type="FORK")
registerDoParallel(cl)

# Kill cluster
on.exit(stopCluster(cl))

# Set a seed
set.seed(1000)


# Function Library --------------------------------------------------------
cluster_hmm <- function(to_cluster){
  cluster_num <- array(NA,length(to_cluster))
  cluster_num [1] <- 1
  for (pp in 2:length (to_cluster)){
    
    if (to_cluster[pp] != to_cluster[(pp-1)]){
      cluster_num[pp] <- cluster_num[(pp-1)] + 1	
    } else {
      cluster_num[pp] <- cluster_num[(pp-1)]	
    }
  }
  cluster_num
}


############################################################################
# # Set dummy variables
# res_dir <- "Helianthus_argophyllus_Todesco_Individual"
# n_cores <- 4
# window_sizes <- c(10,50,100)
# n_states <- c(3,5)

############################################################################
# List all the GEA results files...

gea_res <- list.files(paste0(res_dir,"/"),pattern = "_GEA.rds")
climate_vars <- gsub("_GEA.rds","",gea_res)

# We want to loop over each set off GEA results and output HMM results for each...
for(climate_var in climate_vars){
  
  message(paste0(">>> Starting HMMs for ",climate_var,"
                 "))
  
  # Read in the tau corr
  gea_tau <- readRDS(paste0(res_dir,"/",climate_var,"_GEA.rds"))
  gea_tau <- gea_tau %>% separate(col = "snp_id",into=c("chr","bp"),sep = ":")
  gea_tau$bp <- as.integer(gea_tau$bp)
  
  # Fetch chroms and sizes..
  chr_sizes <- data.frame(gea_tau %>% group_by(chr) %>% summarise(size=max(bp)))
  
  # Cap chromosomes at a megabase...
  chr_sizes <- chr_sizes[chr_sizes$size > 1e6,]
  
  # Save window results
  # wind_results <- array (NA, c((length(window_sizes)*length(which_chroms)*ncol(rho2)),51))
  
  # Loop over chroms
  chr_res <- foreach(chr=chr_sizes$chr) %dopar% {
    message(paste0(">>> Starting ",chr,"
                   "))
    
    # Subset for chrom and take positions
    chrom_r2 <- data.frame(r2=gea_tau[gea_tau$chr == chr,"tau.corr"]^2,
                           bp=gea_tau[gea_tau$chr == chr,"bp"])
    
    # Build windows...
    max_values <- sapply(window_sizes,function(wind) floor(nrow(chrom_r2)/wind))
    for(wind in window_sizes){
      tmp <- rep(seq(1,(nrow(chrom_r2)/wind),1),each=wind)
      chrom_r2$wind <- c(tmp,rep(NA,nrow(chrom_r2)-length(tmp)))
      colnames(chrom_r2)[which(colnames(chrom_r2)=="wind")] <- paste0("wind_",wind)
    }
    
    # Build list of mean rho for each window size
    window_means <- lapply(paste0("wind_",window_sizes),function(wind){
      data.frame(chrom_r2 %>% group_by_at(wind) %>% summarise(mean_r2=mean(r2,na.rm = T),
                                                              var_r2=var(r2,na.rm=T),
                                                              pos=mean(bp)))
    })
    
    # Now we can loop over all window sets and build models...
    message(paste0(">>> Building HMM and Clustering..."))
    window_res <- lapply(window_means,function(input){
      state_res <- lapply(n_states,function(states){
        
        # Fit the model
        mod <- depmix(input[,"mean_r2"] ~ 1,nstates = states, family = gaussian(),ntimes=nrow(input))
        
        # Get fit
        fit.mod <- fit(mod)
        
        # Estimate state
        est.states <- posterior(fit.mod,type="viterbi")
        
        # Cluster states
        est.states$cluster <- cluster_hmm(est.states$state)
        
        # Save the full output
        return(est.states)
      })
      names(state_res) <- paste0(n_states,"_state")
      return(state_res)
    })
    
    # Set names here
    names(window_res) <- paste0("wind_",window_sizes)
    return(window_res)
  }
  
  names(chr_res) <- chr_sizes$chr

# Save these results...
saveRDS(chr_res,paste0(res_dir,"/",climate_var,"_HMM_peak_clusters.rds"))
}
