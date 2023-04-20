# Script for merging gene labels from gff with OF2 IDs
lib <- c("data.table","tidyr","ape")
sapply(lib,library,character.only=T)


# List out all of the gffs and proteomes here... --------------------------
gff_paths <- c("data/reference_genomes/A.halleri_gemmifera_reference/Arabidopsis_halleri.Ahal2.2.45_CHR_MATCHED.gff3",
               "data/reference_genomes/A.lyrata_reference/Alyr_v.1.0_genomic.gff",
               "data/reference_genomes/A.thaliana_reference/GCF_000001735.4_TAIR10.1_genomic_CHR_MATCHED.gff",
               "data/reference_genomes/A.tuberculatus_reference/Amaranthus_tuberculatus_annos1-cds1-id_typename-nu1-upa1-add_chr0.gid54057_CHR_MATCHED.gff",
               "data/reference_genomes/B.stricta_reference/LTM_v2.2_annotatioin_CHR_MATCHED.gff3",
               "data/reference_genomes/C.rubella_reference/GCF_000375325.1_Caprub1_0_genomic_CHR_MATCHED_v2.gff",
               "data/reference_genomes/E.grandis_reference/Egrandis1_0_genomic_CHR_MATCHED.gff",
               "data/reference_genomes/H.annuus_reference/HAN412_Eugene_curated_v1_1.gff3",
               "data/reference_genomes/M.truncatula_reference/Mt4.0v2_genes_20140818_1100_CHR_MATCHED.gff3",
               "data/reference_genomes/P.abies_reference/Pabies01-gene_CHR_MATCHED.gff3",
               "data/reference_genomes/P.deltoides_reference/PdeltoidesWV94_445_v2.1.gene_CHR_MATCHED.gff3",
               "data/reference_genomes/P.hallii_reference/GCF_002211085.1_PHallii_v3.1_genomic_CHR_MATCHED.gff",
               "data/reference_genomes/P.menziesii_reference/Psme.1_0.gtf",
               "data/reference_genomes/P.taeda_reference/Pita.2_01.gtf",
               "data/reference_genomes/P.tremula_reference/Potra02_genes.CHR_MATCHED.gff",
               "data/reference_genomes/P.trichocarpa_reference/Ptrichocarpa_210_v3.0_CHR_MATCHED.gff3",
               "data/reference_genomes/Q.petraea_reference/Qrob_PM1N_genes_20161004.gff",
               "data/reference_genomes/A.alpina_reference/Arabis_alpina.MPIPZ.V5.chr.all.liftOverV4.v3.gff3")

prot_paths <- c("data/proteomes/Ahalleri.faa",
                "data/proteomes/Alyrata.faa",
                "data/proteomes/Athaliana.faa",
                "data/proteomes/Atubercatus.faa",
                "data/proteomes/Bstricta.faa",
                "data/proteomes/Crubella.faa",
                "data/proteomes/Egrandis.faa",
                "data/proteomes/Hannuus.faa",
                "data/proteomes/Mtruncatula.faa",
                "data/proteomes/Pabies.faa",
                "data/proteomes/Pdeltoides.faa",
                "data/proteomes/Phallii.faa",
                "data/proteomes/Pmenziesii.faa",
                "data/proteomes/Ptaeda.faa",
                "data/proteomes/Ptremula.faa",
                "data/proteomes/Ptrichocarpa.faa",
                "data/proteomes/Qpetraea.faa",
                "data/proteomes/Aalpina.faa")

to_transform <- data.frame(gff_path=gff_paths,
                           prot_path=prot_paths)

# Lopp over
for(i in 1:nrow(to_transform)){
  
  # Test subjects
  gff_path <- to_transform$gff_path[i]
  prot_path <- to_transform$prot_path[i]
  
  print(paste0(">>> STARTING ",basename(prot_path)))
  
  # Adjust proteome for merging ---------------------------------------------
  # Fetch the proteome indices
  prot_fai <- system(paste0("grep '>' ",prot_path," | sed 's/>//'"),intern = T)
  
  # Quick check format and kill
  if(length(grep("gene=",prot_fai)) == length(prot_fai)){
    
    print("Already processed, but uncomment if needed...")
    # Separate out the geneID
    prot_fai_split <- strsplit(prot_fai," ")
    prot_gene_ids <- sapply(prot_fai_split,'[[',grep("gene=",prot_fai_split[[1]]))
    prot_gene_dd <- data.frame(gene_ID=gsub("gene=","",prot_gene_ids),
                               OF_ID=(1:length(prot_gene_ids))-1) ## Don't forget there is an out-by-one applied for OF2!!
    
    # Go through and fix any nbis cases...
    to_fix <- grep("nbis",prot_gene_dd$gene_ID)
    prot_gene_dd$gene_ID[to_fix] <- sapply(prot_fai_split[to_fix],'[[',1)
    
    # Adjust gff for merging ---------------------------------------------
    # Set gene boundaries from gff
    if(grepl("gff3",gff_path)){
      #gff <- read.table(gff_path,fill=TRUE,comment.char="#",sep="\t")
      gff <- read.gff(gff_path,GFF3 = TRUE)
    } else {
      gff <- read.gff(gff_path,GFF3 = FALSE)
    }
    #colnames(gff) <- c("seqid","source","type","start","end","score","strand","phase","att")
    
    # Separate out the geneID
    gff_sub <- gff[grep("Parent=",gff$attributes),]
    gff_sub_split_attributes <- strsplit(gff_sub$attributes,";")
    gff_sub$gene_ID <- sapply(gff_sub_split_attributes,function(x) grep("Parent=",x,value=T))
    gff_sub$gene_ID <- gsub("Parent=","",gff_sub$gene_ID)
    gff_sub <- gff_sub[,c(colnames(gff)[1],"start","end","gene_ID")]
    
    # Clean gff to remove duplicated geneIDs and take the largest (which seems to correspond to the "gene" identifier...)
    gff_sub$size <- gff_sub$end - gff_sub$start
    gff_sub <- gff_sub[order(-gff_sub$size),]
    gff_sub_clean <- gff_sub[!(duplicated(gff_sub$gene_ID)),]
    
    # Merge the two together... -----------------------------------------------
    # First check if any missing as a sanity check and just pull the first fai element if they are...
    still_missing <- which(!(prot_gene_dd$gene_ID %in% gff_sub_clean$gene_ID))
    # Finally, grep through any cases of nbis-mrna which doesn't appear at all in the gff...
    if(length(still_missing) > 0){
      for(nbis in still_missing){
        # Get the gene name
        name_check <- gsub("name=","",prot_fai_split[[nbis]][grep("name=",prot_fai_split[[nbis]])])
        chr_check <- gsub("seq_id=","",prot_fai_split[[nbis]][grep("seq_id=",prot_fai_split[[nbis]])])
        gff_tmp <- gff[grepl(name_check,gff$att) & gff[,1] == chr_check,]
        new_name <- gff_tmp[grep("Parent=gene",gff_tmp$attributes),]
        new_name <- unique(sapply(strsplit(new_name$attributes,";"),function(x) grep("Parent=",x,value=T)))
        if(length(new_name)==1){
          print(paste0("Successfully renaming ",prot_gene_dd$gene_ID[nbis]))
          prot_gene_dd$gene_ID[nbis] <- gsub("Parent=","",new_name)
        } else {
          print(paste0("Can't rename ",prot_gene_dd$gene_ID[nbis]))
        }
      }
    }
    
    # And one final check
    still_missing <- which(!(prot_gene_dd$gene_ID %in% gff_sub_clean$gene_ID))
    prot_gene_dd$gene_ID[still_missing] <- sapply(prot_fai_split[still_missing],'[[',1)
    
    # NOW MERGE
    prot_gene_dd$gene_ID <- as.character(prot_gene_dd$gene_ID)
    prot_gff_merge <- merge(prot_gene_dd,gff_sub_clean[,c(colnames(gff)[1],"start","end","gene_ID")],by="gene_ID")
    prot_gff_merge$gea_gene <- paste0(prot_gff_merge[,colnames(gff)[1]],":",prot_gff_merge$start,"-",prot_gff_merge$end)
    
    # Check our output is complete and right to results...
    prop_proteome_matched <- nrow(prot_gff_merge)/length(prot_fai)
    print(paste0(round(prop_proteome_matched*100,3),"% of proteome matched to gff..."))
    
    
    if(prop_proteome_matched==1){
      # Save this to the reference genome directory...
      output_dir <- dirname(gff_path)
      write.table(prot_gff_merge,
                  paste0(output_dir,"/",gsub(".faa","",basename(prot_path)),"_proteome_to_OF_id_map.txt"),
                  row.names = F,quote = F,sep="\t")
    }
    
  } else {
    # OTHERWISE, WE ASSUME THAT THE PROTEOME WAS TAKEN EXTERNALLY BUT IS MATCHED TO THE GFF, SO CAN BE ASSIGNED AS SO...
    print("Proteome/GFF are presumed to have been derived externally...")
    
    # Split protein up and take first element for gene name
    # Run as before mind...
    prot_gene_dd <- data.frame(gene_ID=sapply(strsplit(prot_fai," "),'[[',1),
                               OF_ID=(1:length(prot_fai))-1)
    
    # Resolve the Pabies naming convention...
    if(grepl("Pabies",prot_path)){
      prot_gene_dd$gene_ID <- gsub("MA","Pabies",prot_gene_dd$gene_ID)
    }
    
    # Resolve the Ptremula transcripts issue...
    if(grepl("Ptremula",prot_path)){
      prot_gene_dd$gene_ID <- sapply(strsplit(prot_gene_dd$gene_ID,"\\."),'[[',1)
    }
    
    # Fetch gff
    if(grepl("gff3",gff_path)){
      #gff <- read.table(gff_path,fill=TRUE,comment.char="#",sep="\t")
      gff <- read.gff(gff_path,GFF3 = TRUE)
    } else {
      gff <- read.gff(gff_path,GFF3 = FALSE)
    }
    
    # Set up as genes
    if(grepl("Q.pet",gff_path)){
      gff_sub <- gff[gff[,3]=="mRNA",]
    } else {
      gff_sub <- gff[gff[,3]=="gene",]
    }
    
    if(any(grepl("ID=",gff_sub$attributes))){
      split_attributes <- strsplit(gff_sub$attributes,";")
      gff_sub$gene_ID <- gsub("ID=","",sapply(split_attributes, function(x) grep("ID=",x,value=T)))
    } else {
      gff_sub$gene_ID <- gff_sub$attributes
    }

    # Change transcript to protein to match...
    if(grepl("Q.pet",gff_path)){
      gff_sub$gene_ID <- gsub("T0","P0",gff_sub$gene_ID)
    }
    
    # Merge
    prot_gff_merge <- merge(prot_gene_dd,gff_sub[,c(colnames(gff)[1],"start","end","gene_ID")],by="gene_ID")
    prot_gff_merge$gea_gene <- paste0(prot_gff_merge$seqname,":",prot_gff_merge$start,"-",prot_gff_merge$end)
    
    # Check our output is complete and right to results...
    prop_proteome_matched <- nrow(prot_gff_merge)/length(prot_fai)
    print(paste0(round(prop_proteome_matched*100,3),"% of proteome matched to gff..."))
    
    if(prop_proteome_matched==1){
      # Save this to the reference genome directory...
      output_dir <- dirname(gff_path)
      write.table(prot_gff_merge,
                  paste0(output_dir,"/",gsub(".faa","",basename(prot_path)),"_proteome_to_OF_id_map.txt"),
                  row.names = F,quote = F,sep="\t")
    }
    
  }
  
}


