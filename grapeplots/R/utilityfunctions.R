# Contact: Jimmy  mathewj2@mskcc.org
#
#' @importFrom stats var
#' @import NMI
library(CePa)
library(NMI)
library(mclust)
library(igraph)

#####################
### Function list ###
#####################
#
#  get_sample_data
#  cleaned_gene_names
#  cleaned_sample_names
#  check_numeric
#  create_igraph
#  get_common_names
#  normalize
#  mean_normalize
#  zscore_normalize
#  variance_normalize
#  linkage_index_to_name

get_sample_data <- function(node_data_file) {
    if(grepl("\\.(gct|GCT)$", node_data_file)) {
        data <- t(read.gct(node_data_file))
        colnames(data) <- cleaned_gene_names(colnames(data))
        rownames(data) <- cleaned_sample_names(rownames(data))
        return(data)
    } else {
        data <- read.table(node_data_file, header=TRUE, stringsAsFactors=FALSE, sep=',')
        sample_names <- data[,1]
        tally <- 0
        for(i in 1:length(sample_names)) {
            if(check_numeric(sample_names[i])) {
                tally <- tally + 1
            }
        }
        if(tally > length(sample_names)/2) {
            stop(paste("Check integrity of file", node_data_file, "; majority of sample names are actually numeric."))
        }
        data <- data[,2:dim(data)[2]]
        rownames(data) <- sample_names
        data <- data[,!is.na(colSums(data))]
        data <- data[!is.na(rowSums(data)),]
        return(data)
    }
}

# Special for certain GCT formats
cleaned_gene_names <- function(old_names) {
    for(i in 1:length(old_names)) 
    {
        pattern <- "\\|[0-9]+"
        replacement <- ""
        old_names[i] <- sub(pattern, replacement, old_names[i])
    }
    return(old_names)
}

# Special for certain GCT formats
cleaned_sample_names <- function(old_names) {
    for(i in 1:length(old_names)) 
    {
        pattern <- "\\."
        replacement <- "-"
        old_names[i] <- gsub(pattern, replacement, old_names[i])
    }
    return(old_names)
}

check_numeric <- function(variable) {
    tryCatch(
    {
        if(grepl("^[0-9]+\\.?[0-9]*$", variable)) {
            number_version <- as.numeric(variable)
            if(is.na(number_version)) {
                return(FALSE)
            }
            return(TRUE)
        }
    },
    error=function(err){
        return(FALSE)
    })
    return(FALSE)
}

create_igraph <- function(edge_list_file) {
    if(grepl("\\.(graphml|GRAPHML)$", edge_list_file)) {
        cat("TRY\n")
        graph <- igraph::simplify(read_graph(edge_list_file, format = "graphml"))
        cat("PASS\n")
    } else if(grepl("\\.(ncol|NCOL)$", edge_list_file)) {
        graph <- igraph::simplify(read_graph(edge_list_file, format = "ncol", directed = FALSE))
    } else {
        stop(paste0("Graph file (", edge_list_file, ") format not supported (only graphml and ncol supported)."))
    }
    return(graph)
}

get_common_names <- function(data, graph) {
    node_names <-V(graph)$name
    sample_attribute_names <- colnames(data)

    d1 <- setdiff(node_names, sample_attribute_names)
    d2 <- setdiff(sample_attribute_names, node_names)

    nodes <- intersect(node_names, sample_attribute_names)
    if(length(nodes) <= 3) {
        stop("3 or fewer nodes in the network with data.")
    }
    cat(paste0("      Number of nodes in network not in the data set: ", length(d1), "\n"))
    if(length(d1) > 0) {
        cat(paste0("      ", paste0(d1[1:min(length(d1), 5)], collapse=" ") ))
        if(length(d1) > 5) {
            cat(" ...")
        }
        cat("\n")
    }
    cat(paste0("      Number of nodes in the data set not in the network: ", length(d2), "\n"))
    if(length(d2) > 0) {
        cat(paste0("      ", paste0(d2[1:min(length(d2), 5)], collapse=" ") ))
        if(length(d2) > 5) {
            cat(" ...")
        }
        cat("\n")
    }

    cat(paste0("      Number of nodes in common: ", length(nodes), "\n"))
    return(nodes)
}

normalize <- function(data, normalization_method=c("dividebymean", "zscore", "none")) {
    normalization_method <- match.arg(normalization_method)
    if(normalization_method == "dividebymean") {
        return(mean_normalize(data))
    }
    if(normalization_method == "zscore") {
        return(zscore_normalize(data))
    }
    if(normalization_method == "none") {
        return(data)
    }
    return(NA)   
}

mean_normalize <- function(data) {
    if(!all(data >= 0) || any(is.na(data))) {
        warning("Using divide-by-mean normalization but not all data values are non-negative.")
    }
    for(i in 1:ncol(data)) {
        mean_i <- mean(data[,i])
        data[,i] <- data[,i] / mean_i
    }
    return(data)
}

zscore_normalize <- function(data) {
    if(any(is.na(data))) {
        warning("Some values are NA!")
    }
    for(i in 1:ncol(data)) {
        variance_i <- var(data[,i])
        mean_i <- mean(data[,i])
        if(is.na(variance_i) || variance_i == 0) {
            warning(paste0("Variable ", i, " is constant."))
            data[,i] <- 0
        } else {
            data[,i] <- (data[,i]-mean_i) / sqrt(variance_i)
        }
    }
    return(data)
}

variance_normalize <- function(data) {
    if(!all(data >= 0) || any(is.na(data))) {
        warning("Using divide-by-variance normalization but not all data values are non-negative.")
    }
    for(i in 1:ncol(data)) {
        variance_i <- var(data[,i])
        data[,i] <- data[,i] / sqrt(variance_i)
    }
    return(data)
}

# The linkage function uses a neat but strange indexing system for groups of nodes. The system is unwrapped here.
linkage_index_to_name <- function(index, node_names) {
    if(index < 0) {
        return(node_names[-index])
    }
    if(index > 0) {
        return(paste0("SN", index))
    }
    warning(paste0("Attempted to lookup a node name which does not exist: index ", index, " out of ", length(node_names)))
    return("No name")
}

#' Calculates statistics relating a hierarchical clustering to a known partition
#'
#' Helper function for calling \code{calculate_series_of_nmis_from_analyzed_graph} or \code{calculate_series_of_aris_from_analyzed_graph}.
#' @param comparison_cluster_file A CSV file (no header) with a column of names and a column of partition labels
#' @param graphml_file The GraphML file containing the hierarchy (e.g. as output by \code{generate_reduction}).
#' @param track_partitions (optional) Whether to return a data frame describing the hierarchy of partitions in \code{linkage} (default FALSE)
#' @return NMI and ARI values (possibly also a partitions data frame)
#' @export
calculate_stats_from_graphml_hierarchy <- function(comparison_cluster_file, graphml_file, track_partitions=FALSE) {
    graph <- read_graph(graphml_file, format="graphml")
    return(calculate_series_of_stats_from_analyzed_graph(comparison_cluster_file, graph, track_partitions=track_partitions))
}

#' Calculates statistics relating a linkage/hierarchical clustering to a known partition
#'
#' Helper function for calling \code{calculate_series_of_nmis_from_linkage}.
#' @param comparison_cluster_file A CSV file (no header) with a column of names and a column of partition labels
#' @param graph The graph as returned by \code{hclust}
#' @param track_partitions (optional) Whether to return a data frame describing the hierarchy of partitions in \code{linkage} (default FALSE)
#' @return nmi and ari values (possibly also a partitions data frame)
#' @export
calculate_series_of_stats_from_analyzed_graph <- function(comparison_cluster_file, graph, track_partitions=FALSE) {
    vedges <- E(graph)[E(graph)$virtual_type == "virtual"]
    vgraph <- subgraph.edges(graph, vedges, delete.vertices=FALSE)
    weights <- as_adjacency_matrix(vgraph, attr = "average_distance", sparse = TRUE)
    weights[weights == 0] <- 5*max(weights) # 0 weight only means that there were no relevant triangles # Maybe consider NA here instead?
    # for(i in 1:dim(weights)[1]) {
    #     cat(paste0(i, " nodes.\r"))
    #     weights[i,i] <- 0
    # }
    dist <- as.dist(weights)
    linkage_average <- hclust(dist, method = "average", members = NULL)
    return(calculate_series_of_stats_from_linkage(comparison_cluster_file, linkage_average, track_partitions=track_partitions))
}

#' Calculates statistics relating a linkage/hierarchical clustering to a known partition
#'
#' @param comparison_cluster_file A CSV file (no header) with a column of names and a column of partition labels
#' @param linkage A linkage, as returned by \code{hclust}
#' @param track_partitions (optional) Whether to return a data frame describing the hierarchy of partitions in \code{linkage} (default FALSE)
#' @return nmi values (possibly also a partitions data frame)
#' @export
calculate_series_of_stats_from_linkage <- function(comparison_cluster_file, linkage, track_partitions=FALSE) {
    cf <- read.csv(comparison_cluster_file, stringsAsFactors=FALSE, header=TRUE)
    sample_names <- sapply(cf[,1], function(x){as.character(x)})
    cluster_assignments <- cf[,2]
    # if(length(sample_names)!= length(linkage$labels) || !all(sort(sample_names) == sort(linkage$labels))) {
    #     stop("In comparison clustering, different node set.")
    # }

    names(cluster_assignments) <- sample_names
    matched_cluster_assignments <- sapply(linkage$labels, FUN=function(x){if(x %in% sample_names){return(cluster_assignments[as.character(x)])} else{return("NONE")}})
    names(matched_cluster_assignments) <- linkage$labels

    number_nodes <- length(linkage$labels)
    running_classification <- data.frame(group=(-1*c(1:number_nodes)), row.names=linkage$labels)

    merges <- linkage$merge
    heights <- linkage$height

    number_merges <- dim(merges)[1]
    nmis <- c(rep(-1, number_merges))
    aris <- c(rep(-1, number_merges))

    if(track_partitions) {
        partitions <- list()
        k <- 1
    }

    cat("\n")
    for(i in 1:number_merges) {
        cat(paste0("\rCalculating NMI,ARI after merge ", i, "/", number_merges," ... "))
        name1 <- linkage_index_to_name(merges[i, 1], linkage$labels)
        name2 <- linkage_index_to_name(merges[i, 2], linkage$labels)
        running_classification$group[running_classification$group == merges[i, 1]] <- i
        running_classification$group[running_classification$group == merges[i, 2]] <- i

        # Testing out an in-loop improvement step of the partition
        bestpartition <- running_classification$group
        names(bestpartition) <- rownames(running_classification)
        counts <- sort(table(bestpartition), decreasing=TRUE)
        bestparts <- names(which(counts > 0.01*sum(counts)))
        modifiedpartition <- sapply(bestpartition, FUN=function(x){if(x %in% bestparts){return(x)} else{return(0)}})

        # nmis[i] <- NMI(data.frame(c(1:number_nodes),sorted_assignments),
        #                data.frame(c(1:number_nodes),running_classification))
        # aris[i] <- mclust::adjustedRandIndex(as.vector(sorted_assignments), as.vector(running_classification$group))

        # nmis[i] <- NMI(data.frame(c(1:number_nodes), cluster_assignments),
        #                data.frame(c(1:number_nodes), modifiedpartition[names(cluster_assignments)]))
        # aris[i] <- mclust::adjustedRandIndex(as.vector(cluster_assignments), as.vector(modifiedpartition))
        nmis[i] <- NMI(data.frame(c(1:number_nodes), matched_cluster_assignments),
                       data.frame(c(1:number_nodes), modifiedpartition[names(matched_cluster_assignments)]))
        aris[i] <- mclust::adjustedRandIndex(as.vector(matched_cluster_assignments), as.vector(modifiedpartition[names(matched_cluster_assignments)]))

        if(track_partitions) {
            # partitions[[k]] <- running_classification
            partitions[[k]] <- modifiedpartition
            k <- k+1
        }
        cat(paste0(format(nmis[i], digits=3), ",",format(aris[i], digits=3)))
    }
    cat("\n")
    if(track_partitions) {
        partitions <- data.frame(partitions)
        colnames(partitions) <- c()
        return(list(nmis=nmis, aris=aris, partitions=partitions))
    } else {
        return(list(nmis=nmis, aris=aris))
    }
}   


