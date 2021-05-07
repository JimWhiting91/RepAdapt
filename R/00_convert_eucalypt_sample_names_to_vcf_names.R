# Script renames Eucalypts based on VCF sample names...
lib <- c("data.table")
lapply(lib,library,character.only=T)

# Read in
dd <- read.csv("data/eucalypts_samples_vcf_ids.csv",header=F)

# Fetch names in order
sample_names <- dd[,1]

# Remove whitespace
sample_names <- gsub(" ","",sample_names)

# Ordered VCF names
ordered_vcf_names <- NULL
vcf_names <- dd[,2]
for(i in 1:length(sample_names)){
  print(i)
  
  # Fetch names that match
  tmp <- grep(paste0(sample_names[i],"_"),vcf_names,value=T)
  
  # Take the first one, and do not replace it...
  ordered_vcf_names[i] <- tmp[1]

  # Remove from the set
  if(length(tmp) != 0){
  vcf_names <- vcf_names[vcf_names != tmp[1]]
  }
}

# Write these
write.table(data.frame(ordered_vcf_names),
            "outputs/eucalyptus_vcf_names.txt",
            sep="\t",quote = F,row.names = F,col.names = F)
