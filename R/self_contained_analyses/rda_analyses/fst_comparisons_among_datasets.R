lib <- c("ggplot2","data.table")
sapply(lib,library,character.only=T)

# Read in all the fst results
res_dirs <- list.files("outputs/GEA_res")
res_dirs <- na.omit(sapply(res_dirs,function(dir) ifelse(file.exists(paste0("outputs/GEA_res/",dir,"/SNPRelate_pca_fst_results.rds")),dir,NA)))
fst_res <- lapply(res_dirs,function(dir){
  tmp <- readRDS(paste0("outputs/GEA_res/",dir,"/SNPRelate_pca_fst_results.rds"))
})

all_fst <- data.frame(rbindlist(lapply(names(fst_res),function(dataset){
  print(dataset)
  if(!(is.null(fst_res[[dataset]]$fst$MeanFst))){
    data.frame(dataset=dataset,
               fst=fst_res[[dataset]]$fst$MeanFst,
               eig1=fst_res[[dataset]]$eig1,
               eig50=fst_res[[dataset]]$eig50/length(fst_res[[dataset]]$pca$eigenval))
  } else {
    data.frame(dataset=dataset,
               fst=fst_res[[dataset]]$fst$FST,
               eig1=fst_res[[dataset]]$eig1,
               eig50=fst_res[[dataset]]$eig50/length(fst_res[[dataset]]$pca$eigenval))
  }
})))

ggplot(all_fst,aes(x=fst,y=eig50))+
  geom_point()+
  geom_text(aes(label=dataset))
