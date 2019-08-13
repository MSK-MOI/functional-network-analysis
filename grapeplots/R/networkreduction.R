# Contact: Jimmy  mathewj2@mskcc.org

#' @import igraph
#' @import parallel
#' @import pbmcapply
#' @import CePa
#' @importFrom stats as.dist hclust complete.cases na.omit cor
#' @importFrom utils read.table read.csv

suppressMessages(library(igraph))
library(parallel)
library(pbmcapply)
library(stats)
library(utils)
library(CePa)

#' Network reduction hierarchy generator
#' 
#' @param node_data_file Data file with numerical node sample data, with node names. Can be simple CSV, or GCT format.
#' @param topology_file Network data file as edge list in \code{.ncol} format, 'nodename1 nodename2 weight', or else \code{.graphml} file format.
#' @param correlation_cutoff (optional) Pearson correlation cutoff for inferring network topology from sample set. Default 0.75.
#' @param connectivity_cutoff (optional) Edge connectivity cutoff (also based on Pearson correlation) for inferring network topology from sample set.
#' @param iterations (optional) Number of iterations for graph-Laplacian-based node data synthesis (when node_data_file is not specified; default 50 in this case).
#' @param edge_limit (optional) A hard limit on the number of edges in the resulting graph (using any method).
#' @param number_of_gmm_populations (optional) Number of populations for Gaussian Mixture Modeling (default 3).
#' @param correlation_transformation (optional) Transform feature vectors into vector of Pearson correlations against other features. If selected, no normalization is performed beforehand.
#' @param mc.cores (optional) For parallelization, requested number of cores.
#' @param normalization (optional) One of "mean", "zscore", "variance", or "none". Default "none".
#' @param method (optional) Default "gmt" (Gaussian Mixture Transport). If "correlation", skips the whole algorithm and simply calculates correlation weights on original network to create hierarchy.
#' @param logdir (optional) Directory for verbose output log file to be stored in.  Default NA, with output to standard out.
#' @return named list containing: \code{graph} (with weighted virtual edges); \code{hierarchy}, the tree graph encoding the hierarchical clustering of features/nodes with respect to the weights as distance; \code{virtual_edges_and_weights}, a data frame listing the virtual edges (neighbor-of-neighbor) with average GMT distance weights; \code{descriptor}, a string describing most of the parameters of the run. Also writes \code{graph} and \code{hierarchy} to 2 GraphML files, with names derived from \code{descriptor}.
#' @export
generate_reduction <- function(node_data_file=NA,
                               topology_file=NA,
                               correlation_cutoff=0.75,
                               connectivity_cutoff=NA,
                               iterations=NA,
                               edge_limit=NA,
                               number_of_gmm_populations=3,
                               correlation_transformation=c("none", "pearson", "spearman"),
                               normalization=c("none", "zscore", "dividebymean"),
                               method=c("gmt", "simplecorrelation"),
                               mc.cores=1,
                               logdir=NA
                               ) {
    correlation_transformation <- match.arg(correlation_transformation)
    normalization <- match.arg(normalization)
    method <- match.arg(method)

    descriptor <- get_file_descriptor(node_data_file=node_data_file,
                                      topology_file=topology_file,
                                      iterations=iterations,
                                      correlation_transformation=correlation_transformation,
                                      normalization=normalization,
                                      method=method)

    if(is.na(logdir)) {
        log_place <- stdout()
        named_file <- FALSE
    } else {
        log_place <- paste(logdir, paste0(descriptor, ".txt"), sep="/")
        named_file <- TRUE
    }

    cat("\n", file=log_place, append=named_file)
    cat(paste0("Calculating weighted network reduction based on\n",
        "\nnode_data_file:             ", node_data_file,
        "\ntopology_file:              ", topology_file,
        "\ncorrelation_transformation: ", correlation_transformation,
        "\nnormalization:              ", normalization,
        "\nmethod:                     ", method,
        "\n"), file=log_place, append=named_file)
    start_time <- Sys.time()

    cat("\n", file=log_place, append=named_file)
    cat("[1/6] Loading data and loading/building network\n", file=log_place, append=named_file)
    if(is.na(node_data_file) && is.na(topology_file)) {
        stop("      Must specify either a node_data_file or topology_file.")
    }

    if(!is.na(node_data_file)) { data <- get_sample_data(node_data_file) }
    if(!is.na(topology_file)) { graph <- create_igraph(topology_file) }

    if(is.na(node_data_file)) {
        cat(paste0("      Synthesizing numerical data along nodes of network. (No node_data_file supplied).\n"), file=log_place, append=named_file)
        if(is.na(iterations)) {
            iterations <- 50
        }
        data <- synthesize_node_data(graph, iterations=iterations)
    }

    if(correlation_transformation != "none") {
        data <- cor(data, method=correlation_transformation)
    } else {
        data <- normalize(data, normalization_method=normalization)
    }

    if(is.na(topology_file)) {
        cat(paste0("      Inferring feature network from data using Pearson correlation. (No topology_file supplied).\n"), file=log_place, append=named_file)
        graph <- infer_igraph(data, correlation_cutoff, connectivity_cutoff, edge_limit)        
    }

    cat("\n", file=log_place, append=named_file)
    cat("[2/6] Integrating sample set and network data\n", file=log_place, append=named_file)
    common_names <- get_common_names(data, graph)
    data <- data[, common_names]
    graph <- induced_subgraph(graph, common_names)
    gc()
    edge_list <- as.list(data.frame(t(ends(graph, E(graph))), stringsAsFactors=FALSE))
    cat(paste0("      ", length(edge_list), " edges after integrating dataset and feature network.\n"), file=log_place, append=named_file)


    if(method == "gmt") {
        cat("\n", file=log_place, append=named_file)
        cat("[3/6] Fitting Gaussian mixture models\n", file=log_place, append=named_file)
        core_number <- parallel::detectCores()
        if(is.na(core_number)) {
            cat(paste0("      Trying ", mc.cores, " cores, as requested.\n"), file=log_place, append=named_file)
        } else {
            cat(paste0("      Number of cores according to parallel::detectCores(): ", core_number, "\n"), file=log_place, append=named_file)
            cat(paste0("      Trying ", max(mc.cores, core_number-1), ".\n"), file=log_place, append=named_file)
            mc.cores <- max(mc.cores, core_number-1)
        }

        edge_models <- pbmclapply(edge_list,
                                    function(edge) {
                                        gmm_model_an_edge(edge, data, number_of_gmm_populations = number_of_gmm_populations)
                                    }, mc.cores = mc.cores
                                ) # (Uses parallelization. Ordinary one is lapply. 'mc' stands for 'multi-core'.)
        edge_models <- edge_models[!is.na(edge_models) & !is.null(edge_models)]
        cat(paste0("      ", length(edge_models), " 2D models with ", number_of_gmm_populations, " populations.\n"), file=log_place, append=named_file)


        cat("\n", file=log_place, append=named_file)
        cat("[4/6] Comparing all adjacent models\n", file=log_place, append=named_file)
        edge_pair_statistics <- calculate_edge_pair_distances(graph, edge_models, number_of_gmm_populations = number_of_gmm_populations, mc.cores=mc.cores)
        edge_pair_statistics <- edge_pair_statistics[complete.cases(edge_pair_statistics),]
        rownames(edge_pair_statistics) <- c()
        cat(paste0("      ", dim(edge_pair_statistics)[1], " adjacent edge pairs considered.\n"), file=log_place, append=named_file)


        cat("\n", file=log_place, append=named_file)
        cat("[5/6] Averaging over intermediating triangles to get virtual edges with weights\n", file=log_place, append=named_file)
        virtual_edges_and_weights <- calculate_average_distance_over_triangles(edge_pair_statistics)
        vgraph <- graph_from_data_frame(virtual_edges_and_weights, directed=FALSE, vertices = as.data.frame(common_names))
        cat(paste0("      ", dim(virtual_edges_and_weights)[1], " total virtual edges.\n"), file=log_place, append=named_file)


        cat("\n", file=log_place, append=named_file)
    #    weights <- as_adjacency_matrix(vgraph, attr = "average_distance", sparse = TRUE) #sparse=TRUE needed for large cases?
        E(vgraph)$average_distance <- (1.5/max(E(vgraph)$average_distance)) * E(vgraph)$average_distance
        weights <- as_adjacency_matrix(vgraph, attr = "average_distance", sparse = FALSE)
        weights[weights == 0] <- 1.5*max(weights) # 0 weight only means that there were no relevant triangles # Maybe consider NA here instead?
        for(i in 1:dim(weights)[1]) {
            weights[i,i] <- 0
        }
        dist <- as.dist(weights)
    } else if(method == "simplecorrelation") {
        cat("[3/6] (skipped)\n", file=log_place, append=named_file)
        cat("[4/6] (skipped)\n", file=log_place, append=named_file)
        cat("[5/6] (skipped)\n", file=log_place, append=named_file)
        cat("[6/6] Calculating hierarchy/linkage\n", file=log_place, append=named_file)
        dist <- as.dist(cor(data, method="pearson"))
    }
    
    E(graph)$original_type <- c("original")
    E(vgraph)$virtual_type <- c("virtual")
    fullgraph <- graph + vgraph
    E(fullgraph)$original_type[is.na(E(fullgraph)$original_type)] <- c("not original")
    E(fullgraph)$virtual_type[is.na(E(fullgraph)$virtual_type)] <- c("not virtual")
    E(fullgraph)$type <- paste(E(fullgraph)$original_type, E(fullgraph)$virtual_type, sep='/')
    E(fullgraph)$average_distance[is.na(E(fullgraph)$average_distance)] <- c(-1)

    linkage_average <- hclust(dist, method = "average", members = NULL)
    hierarchy <- make_graph_from_linkage(linkage_average, common_names, label="linkage_average")

    results <- list(graph=fullgraph,
                    hierarchy=hierarchy,
                    virtual_edges_and_weights=virtual_edges_and_weights)
    cat(paste0("DONE\n"), file=log_place, append=named_file)
    end_time <- Sys.time()
    cat(end_time-start_time, file=log_place, append=named_file)
    cat("\n\n", file=log_place, append=named_file)
    return(results)
}

#' Get descriptor string for a network reduction run
#' 
#' @param node_data_file Data file with numerical node sample data, with node names. Can be simple CSV, or GCT format.
#' @param topology_file Network data file as edge list in \code{.ncol} format, 'nodename1 nodename2 weight', or else \code{.graphml} file format.
#' @param correlation_transformation (optional) Transform feature vectors into vector of Pearson correlations against other features. If selected, no normalization is performed beforehand.
#' @param iterations (optional) Number of iterations for graph-Laplacian-based node data synthesis.
#' @param normalization (optional) One of "mean", "zscore", "variance", "none". Default "mean" (anticipating positive data, to preserve positivity).
#' @param method (optional) Default "gmt" (Gaussian Mixture Transport). If "correlation", skips the whole algorithm and simply calculates correlation weights on original network to create hierarchy.
#' @return String that will be used as the descriptor (file name) of the output of \link{generate_reduction}.
#' @export
get_file_descriptor <- function(
        node_data_file,
        topology_file,
        iterations=NA,
        correlation_transformation=c("none", "pearson", "spearman"),
        normalization=c("dividebymean", "zscore", "none"),
        method=c("gmt", "simplecorrelation") ) {
    correlation_transformation <- match.arg(correlation_transformation)
    normalization <- match.arg(normalization)
    method <- match.arg(method)

    if(is.na(topology_file)) {
        topology_descriptor <- "corrinferredtopology"
    } else {
        topology_descriptor <- sub("\\.(ncol|NCOL|graphml|GRAPHML)$", "", basename(topology_file))
    }

    if(correlation_transformation != "none") {
        preprocessing_descriptor <- "corrtrans"
    } else {
        if(normalization == "none") {
            preprocessing_descriptor <- "unnormalized"
        } else {
            preprocessing_descriptor <- normalization
        }
    }

    if(is.na(node_data_file)) {
        data_descriptor <- "synthesizeddata"
    } else {
        data_descriptor <- sub("\\.(gct|GCT|csv|CSV)$", "", basename(node_data_file))
    }

    if(is.na(iterations)) {
        iteration_descriptor <- ""
    } else {
        iteration_descriptor <- paste0("_", iterations, "iterations")
    }

    return(paste0(data_descriptor, "_", topology_descriptor, "_", method, "_", preprocessing_descriptor, iteration_descriptor))
}





