# Plot a random graph and highlight nodes with high centrality metrics...
lib = c("ggplot2",'igraph','data.table','ggnetwork','network','sna','tidygraph','ggraph')
sapply(lib,library,character.only = T)

# Make a random network with 100 nodes
set.seed(1000)
network_obj <- watts.strogatz.game(1, 20, 2, 0.2, loops = FALSE, multiple = FALSE)

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

