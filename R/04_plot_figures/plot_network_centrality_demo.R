# # Build a matrix to plot
# network_matrix = matrix(ncol = 10, nrow = 10)
# rownames(network_matrix) <- colnames(network_matrix) <- paste0("G",1:10)
# 
# # Fill each row
# network_matrix[1,] = c(0,6,3,1,0,0,0,0,0,0)
# network_matrix[2,] = c(6,0,3,0,2,0,0,0,0,0)
# network_matrix[3,] = c(3,3,0,2,3,6,0,0,0,0)
# network_matrix[4,] = c(1,0,2,0,0,4,0,0,0,0)
# network_matrix[5,] = c(0,2,3,0,0,0,0,0,0,0)
# network_matrix[6,] = c(0,0,6,4,0,0,2,0,0,0)
# network_matrix[7,] = c(0,0,0,0,0,2,0,8,8,0)
# network_matrix[8,] = c(0,0,0,0,0,0,8,0,8,0)
# network_matrix[9,] = c(0,0,0,0,0,0,8,8,0,0)
# network_matrix[10,] = c(0,0,0,0,0,0,0,0,0,0)
# 
# # # Also make a matrix where only top half...
# # network_matrix2 = network_matrix
# # network_matrix2[lower.tri(network_matrix2)] <- NA
# 
# # Visualise...
# # build the graph object
# library(igraph)
# network_to_plot = data.table(na.omit(melt(network_matrix)))
# network_to_plot = network_to_plot[value > 0,]
# network_obj <- graph_from_edgelist(as.matrix(network_to_plot[,1:2]))
# E(network_obj)$weight = network_to_plot$value
# 
# # plot it
# plot(network_obj)
# 
# # Fetch network stats...
# # Betweenness, cutoff of -1 ensures that no cutoff is used for betweeness
# centrality_stats = data.table(Betweenness = estimate_betweenness(network_obj,cutoff = -1,directed = F),
#                               Closeness = closeness(network_obj,normalized = T),
#                               Strength = sapply(paste0("G",1:9),function(x) sum(network_to_plot[Var1 == x | Var2 == x,value])),
#                               Degree = sapply(paste0("G",1:9),function(x) length(network_to_plot[Var1 == x | Var2 == x,value])),
#                               Node = paste0("G",1:9))
# centrality_fig = melt(centrality_stats) |>
#   ggplot(aes(y = value, x = Node,fill = Node)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(~variable,ncol=2,scale = "free_y",strip.position = 'top') +
#   theme_minimal() + 
#   theme(axis.text.x = element_blank(),
#         axis.title = element_blank(),
#         strip.text = element_text(size = 10),
#         panel.grid.major.x = element_blank(),
#         legend.position = "none")+
#   ggtitle("Centrality Measures") +
#   scale_fill_brewer(palette = "Spectral")
# 
# library(ggnetwork)
# set.seed(10020)
# 
# flattened_network = ggnetwork(network_obj)
# network_fig = ggplot(flattened_network, aes(x = x, y = y, xend = xend, yend = yend)) +
#   geom_edges(aes(alpha = weight),color = "grey50", curvature = 0,size = 2) +
#   theme_blank() +
#   geom_nodelabel(fontface = "bold",aes(label = name, fill = name),size = 5)+
#   theme(legend.position = "none") +
#   scale_fill_brewer(palette = "Spectral")
# 
# 
# 
# centrality_combined = cowplot::plot_grid(centrality_fig,
#                                          network_fig,ncol=2)
# centrality_combined

# Plot a random graph and highlight nodes with high centrality metrics...
lib = c("ggplot2",'igraph','data.table','ggnetwork','network','sna','tidygraph','ggraph')
sapply(lib,library,character.only = T)

# Make a random network with 100 nodes
set.seed(1000)
network_obj <- watts.strogatz.game(1, 20, 2, 0.2, loops = FALSE, multiple = FALSE)
# network_obj = network(rgraph(100, tprob = 0.01), directed = FALSE)
# network_to_plot = data.table(na.omit(melt(network_matrix)))
# network_to_plot = network_to_plot[value > 0,]
# network_obj <- graph_from_edgelist(as.matrix(network_to_plot[,1:2]))
# E(network_obj)$weight = network_to_plot$value

# plot it
plot(network_obj)

# Fetch network stats...
# Betweenness, cutoff of -1 ensures that no cutoff is used for betweeness
centrality_stats = data.frame(Betweenness =  igraph::estimate_betweenness(network_obj,cutoff = -1,directed = F),
                              Closeness = igraph::closeness(network_obj,normalized = T),
                              Degree =  igraph::degree(network_obj))
centrality_stats$Strength = sapply(centrality_stats$Degree,function(x) sum(abs(rnorm(x))))

# Convert
tidy_network <- as_tbl_graph(network_obj) |>
  activate(nodes) |>
  mutate(Betweenness = centrality_stats$Betweenness,
         Closeness = centrality_stats$Closeness,
         Degree = centrality_stats$Degree,
         Strength = centrality_stats$Strength)

# Plot each of these...
centrality_examples = lapply(1:4,function(x){
  
  plot_network <- as_tbl_graph(network_obj) |>
    activate(nodes) |>
    mutate(to_plot = centrality_stats[,x])
  
  ggraph(plot_network, layout = "auto") + 
    geom_edge_link(colour = 'grey75') + 
    geom_node_point(aes(size = to_plot^3,colour = to_plot^3)) + 
    theme_void() + 
    scale_color_gradient(low = "blue4", high = "red2") +
    # scale_colour_viridis(option = 'C') +
    theme(legend.position = 'none',
          title = element_text(family = 'sans',size = 12)) +
    ggtitle(colnames(centrality_stats[x]))
  
})
centrality_combined = cowplot::plot_grid(plotlist = centrality_examples,ncol = 2)
centrality_combined

