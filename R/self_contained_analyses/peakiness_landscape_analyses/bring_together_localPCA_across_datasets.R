################################################################################
# Read in the localPCA results, identify outliers, and assess for each dataset the proportion of genome marked as outlier
################################################################################
lib <- c("Routliers","ggridges","data.table","ggplot2","dplyr","parallel","Rfast","cowplot")
lapply(lib,library,character.only=T)

# Find all the results...
localPCA_res <- list.files("outputs/localPCA_results/",pattern = "_localPCA_results_windsize100.rds")
datasets <- gsub("_localPCA_results_windsize100.rds","",localPCA_res)

# For each, read in the results and chr-start-end and calculate the proportion "outlying"
outlier_prop <- data.frame(rbindlist(pbmcapply::pbmclapply(datasets,function(dataset){
  
  # Fetch res
  tmp_res <- readRDS(paste0("outputs/localPCA_results/",dataset,"_localPCA_results_windsize100.rds"))
  
  # Separate
  start_end <- data.frame(rbindlist(tmp_res[[1]]))
  localpca_tmp <- tmp_res[[2]]
  
  # Loop over all chrs and identify outliers...
  chr_outliers <- sapply(1:length(localpca_tmp),function(i){
    
    chr_mds <- data.frame(mds1=localpca_tmp[[i]]$points[,1],
                          mds2=localpca_tmp[[i]]$points[,2])
    
    # Get outliers
    mds1_outliers <- outliers_mad(chr_mds$mds1)
    mds2_outliers <- outliers_mad(chr_mds$mds2)
    all_outliers <- c(mds1_outliers$outliers_pos,mds2_outliers$outliers_pos)
    
    if(length(all_outliers) > 0){
      length_outlying <- sum(rowSums(tmp_res[[1]][[i]][all_outliers,c("start","end")]))
    } else {
      length_outlying <- 0
    }
    
    return(length_outlying)
    
  })
  
  # Total proportion outlying
  total_outlying <- sum(chr_outliers)/sum(rowSums(start_end[,c("start","end")]))
  return(data.frame(dataset=dataset,
                    prop_genome_outlier=total_outlying))
},mc.cores=4)))

