# Comparison of OrthoFinder and TimeTree phylogenies for reference genome species
lib = c('ggplot2','ggtree','data.table','ape','phytools')
sapply(lib,library,character.only = T)

# Fetch both of the trees
OF_tree = read.tree('outputs/orthology/Results_221213_18_genomes_Ptaeda_isoforms_removed/Species_Tree/SpeciesTree_rooted.txt')
OF_tree = drop.tip(OF_tree,'Pmenziesii')
TT_tree = read.tree('metadata/230221_genome_species.nwk')

# Plot each
plot(OF_tree)
plot(TT_tree)

# Rename the tips in the OF tree to give full names...
OF_tips = OF_tree$tip.label
OF_tips = gsub("Atubercatus","Atuberculatus",OF_tips)
OF_tips = gsub("Qpetraea","Qrobur",OF_tips)
OF_tree$tip.label = sapply(sub('.', '', OF_tips),grep,x = TT_tree$tip.label,value = T)

# Tidy each of these...
OF_tree$tip.label = gsub("_"," ",OF_tree$tip.label)
TT_tree$tip.label = gsub("_"," ",TT_tree$tip.label)

# Plot together
t1 = OF_tree
t2 = TT_tree
obj<-cophylo(t1,t2)

# Save this
pdf("figs/FigureSX_OF_TT_tree_agreement.pdf",height=5,width = 10)
plot(obj,link.lwd=2,link.lty="dotted",
     link.col=make.transparent("red",0.5))
dev.off()
