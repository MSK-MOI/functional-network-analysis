
suppressMessages(library(igraph))
library(fna)

node_data_file <- c("example_data/lung_tissue_expression_gtex_abridged_100.csv")
normalization <- "none"
topology_file <- c(NA)

results <- generate_reduction(node_data_file, topology_file=topology_file, correlation_cutoff=0.65, verbose=TRUE) # log_file="log.txt"
write_graph(results$graph, "example_data/lung_gtex_weighted_100.graphml", format = "graphml")
write_graph(results$hierarchy, "example_data/lung_gtex_hierarchy_100.graphml", format = "graphml")

# Requires download of goa_human.gaf
# gaf_annotate_graphml("example_data/lung_gtex_hierarchy_100.graphml", verbose=TRUE)






