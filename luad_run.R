
suppressMessages(library(igraph))
library(fna)

node_data_file <- c("example_data/tcga_luad.cor_reduced_3000.data_RNA_Seq_v2_mRNA_median_Zscores.csv")
normalization <- "none"
topology_file <- c(NA)

results <- generate_reduction(node_data_file, topology_file=topology_file, correlation_cutoff=0.55, verbose=TRUE)
write_graph(results$graph, "example_data/luad_weighted_3000.graphml", format = "graphml")
write_graph(results$hierarchy, "example_data/luad_hierarchy_3000.graphml", format = "graphml")

node_data_file_for_clustering <- "example_data/LUAD-TP.normalized.gct"
cluster_file <- "example_data/LUAD-TP.bestclus.csv"
annotated <- annotate_with_clusters(
    results$graph,
    node_data_file_for_clustering,
    cluster_file,
    correlation_transformation="none",
    normalization="none",
    verbose=TRUE
)
write_graph(annotated, "example_data/luad_weighted_3000_gdac_annotated.graphml", format = "graphml")

# Requires download of goa_human.gaf
#gaf_annotate_graphml("example_data/lung_gtex_hierarchy_100.graphml", verbose=TRUE)






