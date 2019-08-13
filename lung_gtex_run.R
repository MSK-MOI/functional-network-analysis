
suppressMessages(library(igraph))
library(grapeplots)

node_data_file <- c("example_data/lung_tissue_expression_gtex_abridged.csv")
normalization <- "none"
topology_file <- c(NA)

results <- generate_reduction(node_data_file, topology_file=topology_file, correlation_cutoff=0.7)
write_graph(results$graph, "example_data/lung_gtex_weighted.graphml", format = "graphml")
write_graph(results$hierarchy, "example_data/lung_gtex_hierarchy.graphml", format = "graphml")

graph<-read_graph("example_data/lung_gtex_weighted.graphml", format = "graphml")
hierarchy<-read_graph("example_data/lung_gtex_hierarchy.graphml", format = "graphml")

