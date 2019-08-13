# Contact: Jimmy  mathewj2@mskcc.org

#' Synthesize node data from a graph topology
#'
#' @param graph igraph object with the topology.
#' @param iterations Number of iterations for smoothing operator. Should probably be selected in terms of average degree, connectivity, etc. Default 50.
#' @param seeds Number of 'seed' node locations to choose randomly for each sample. Default is 5 * the number of nodes divided by the number of samples.
#' @param trials Number of samples (each sample is a numerical vector over the node set). Default is 1/2 the number of nodes.
#' @return A data frame, formatted in the same manner as the ouptut of get_sample_data.
#' @export
synthesize_node_data <- function(graph, iterations=50, seeds=NA, trials=NA) {
    number_nodes <- length(V(graph))
    if(is.na(trials)) {
        trials <- as.integer(0.5 * number_nodes)
    }
    if(is.na(seeds)) {
        seeds <- as.integer(5 * number_nodes / trials)
    }

    lap <- laplacian_matrix(graph, normalized=TRUE)
    initial_distributions <- matrix(0, ncol=trials, nrow=number_nodes)
    initial_distributions <- apply(initial_distributions, 2,
        function(col) {
            ss <- sample(number_nodes, seeds, replace=FALSE)
            return(sapply(c(1:number_nodes), function(node) { as.integer(node %in% ss) } ))
        } )

    # Unclear if it would be faster to precompute lap^iterations, since lap is a sparse matrix (see igraph::laplacian_matrix).
    lap_power <- lap
    for(i in 1:iterations) {
        lap_power <- lap_power %*% lap
    }

    samples <- lap_power %*% initial_distributions
    samples <- t(as.matrix(samples))
    rownames(samples) <- paste0("Sample", c(1:trials))
    if("name" %in% list.vertex.attributes(graph)) {
        colnames(samples) <- V(graph)$name
    }
    return(as.data.frame(samples))
}

#' Infer graph topology from data set
#' 
#' Uses a correlation cutoff (or an equivalent edge connectivity cutoff) to make graph with edges corresponding only to nodes with highly correlated data values
#' @param data The node-data data frame
#' @param correlation_cutoff (optional) A Pearson correlation cutoff value, from 0 to 1. Default 0.75
#' @param connectivity_cutoff (optional) A desired edge connectivity for the resulting graph, from 0 to 1.
#' @param edge_limit (optional) A hard limit on the number of edges that will be present in the resulting graph. Default 100000, set to NA if you want no limit
#' @param log_place (optional) A file/connection to log messages. Default \code{stdout()}.
#' @return An igraph object whose node names come from column names in \code{data}
infer_igraph <- function(data, correlation_cutoff=0.75, connectivity_cutoff=NA, edge_limit=100000, log_place=stdout()) {
    cat(paste0("      Data set ", dim(data)[1], "x", dim(data)[2]), file=log_place, append=TRUE)
    if(all(apply(data, MARGIN = c(1,2), FUN = function(x){is.numeric(x) && !is.na(x)}))) {
        cat(paste0(" all numeric.\n"), file=log_place, append=TRUE)
    } else {
        stop("Data set improperly formatted; not all numeric, or with NAs.")
    }

    variances <- sapply(1:dim(data)[2], FUN=function(i){return(var(data[,i]))})
    number_of_trivial_variables <- sum(variances == 0)
    if(number_of_trivial_variables > 0) {
        cat(paste0("      Found ", number_of_trivial_variables, " trivial variables (with zero variance):\n"))
        cat(paste0(c("      ", colnames(data)[variances == 0], "\n"), collapse=" "))
    }
    data <- data[, variances!=0]

    correlations <- abs(cor(data, method="pearson", use="complete.obs"))
    cat(paste0("      Calculated correlations (", dim(correlations)[1], "x", dim(correlations)[2], ")\n"), file=log_place, append=TRUE)
    number_nodes <- dim(correlations)[1]
    if(!is.na(correlation_cutoff) && !is.na(connectivity_cutoff)) {
        stop("Cannot specify both a correlation cutoff and a connectivity cutoff.\n")
    }
    for(i in 1:number_nodes) {
        for(j in i:1) {
            correlations[i,j] <- 0
        }
    }

    cutoff_value <- 0
    if(!is.na(correlation_cutoff)) {
        if(correlation_cutoff < 0 || correlation_cutoff > 1) {
            stop("Need to supply correlation_cutoff in correct range (0,1].")
        }
        cutoff_value <- correlation_cutoff
    } else {
        if(connectivity_cutoff < 0 || connectivity_cutoff > 1) {
            stop("Need to supply connectivity_cutoff in correct range (0,1].")
        }
        values <- sort(as.vector(correlations), reverse<-TRUE)
        cutoff_index <- as.integer((number_nodes-1) * number_nodes * 0.5 * connectivity_cutoff)
        cutoff_value <- values[cutoff_index]
    }

    cat(paste0("      Cutoff value for correlation: ", cutoff_value, "\n"), file=log_place, append=TRUE)
    edge_list <- c(rep("", 2*(number_nodes-1)*number_nodes/2))
    index <- 1
    if(any(is.na(correlations))) {
        warning("Correlation matrix has NA values.")
    }
    adjacency <- (correlations > cutoff_value)
    if(any(is.na(adjacency))) {
        stop("The matrix of comparisons between correlation matrix entries and the cutoff has NA values.")
    }

    start_time <- Sys.time()
    cat("\n")    
    for(i in 1:number_nodes) {
        stop_time <- Sys.time()
        elapsed <- format(round(difftime(stop_time, start_time, 2, units="mins"), 2), nsmall = 2)
        estimate <- format(round((number_nodes/i)*difftime(stop_time, start_time, units="mins"), 2), nsmall=2) 
        cat(paste0(i, " nodes processed of ", number_nodes, "(Elapsed ", elapsed, " of expected total ", estimate, ")\r"))

        for(j in i:number_nodes){
            if(adjacency[i,j] && j>i) {
                edge_list[index] <- colnames(data)[i]
                edge_list[index+1] <- colnames(data)[j]
                index <- index+2
            }
            if(!is.na(edge_limit)) {
                if(index > 2*edge_limit) {
                    break
                }
            }
        }
        if(!is.na(edge_limit)) {
            if(index > 2*edge_limit) {
                cat(paste0("      Reached edge_limit ", edge_limit, ".\n"), file=log_place, append=TRUE)
                break
            }
        }
    }
    cat("\n")
    edge_list <- edge_list[1:(index-1)] # Since the correct number of edges might not be *exactly* as predicted by the calculation of cutoff_index
    n <- length(edge_list)/2
    cat(paste0("      Inferred ", n, " edges, out of possible ", (number_nodes-1) * number_nodes *0.5, ", with connectivity ", n/((number_nodes-1)*number_nodes * 0.5), "\n" ), file=log_place, append=TRUE)
    inferred <- make_graph(edge_list, directed=FALSE)
    cat(paste0("      Support of inferred graph contains ", length(V(inferred)), " nodes.\n"), file=log_place, append=TRUE)
    cat(paste0("      Average degree ", mean(degree(inferred)), ".\n"), file=log_place, append=TRUE)
    cat(paste0("      Diameter ", diameter(inferred), ".\n"), file=log_place, append=TRUE)
    return(inferred)
}

#' Tabulation of W2 distances between GMMs for edges of a graph
#' 
#' @param graph Given igraph graph object
#' @param edge_models List of GMMs
#' @param number_of_gmm_populations (optional) Number of populations for the Gaussian Mixture Modeling (default 3)
#' @param mc.cores (optional) Number of cores requested for parallelization (default 1)
#' @return A data frame with the two edge indices, the distance, and the three endpoint nodes involved
calculate_edge_pair_distances <- function(graph, edge_models, number_of_gmm_populations=3, mc.cores=1) {
    node_index <- as.matrix(V(graph))
    edges <- ends(graph, E(graph))
    index_adjacency <- as_adjacency_matrix(graph, edges = TRUE, sparse = TRUE)

    # Uses parallelization
    edge_pair_statistics <- pbmclapply(c(1:length(node_index)),
        function(k) {
            calculate_edge_pair_distances_based_at(k, index_adjacency, edge_models, edges, number_of_gmm_populations)
        }, mc.cores=mc.cores) # If warnings (or errors?) occur inside the pbmclapply call, sometimes they do not break the program flow, which continues to the below code with junk data. This has something to do with forking processes or threads or something.

    # Unparallelized version
    # edge_pair_statistics <- lapply(c(1:length(node_index)),
    #     function(k) {
    #         calculate_edge_pair_distances_based_at(k, index_adjacency, edge_models, edges, number_of_gmm_populations)
    #     }) 

    eps <- data.frame()
    for(k in 1:length(node_index)) {
        if(length(edge_pair_statistics) == 0 ) {
            stop("Bad error, edge_pair_statistics is empty.")
        }
        cond <- ( length(edge_pair_statistics[[k]]) > 0 )
        if(cond) {
            eps <- rbind(eps, edge_pair_statistics[[k]])
        }
    }
    eps <- na.omit(eps)
    # The following sort order ensures that all "triangles" based on a given pair of nodes are consecutive
    dictionary_order <- order(eps[, "endpoint1"], eps[, "endpoint2"])
    eps <- eps[dictionary_order, ]
    edge_pair_statistics <- eps
    rownames(edge_pair_statistics) <- c()
    return(edge_pair_statistics)
}

calculate_edge_pair_distances_based_at <- function(k, index_adjacency, edge_models, edges, number_of_gmm_populations) {
    stats_k <- list()
    count <- 1
    edges_from_k <- index_adjacency[k, index_adjacency[k, ] != 0] # Gets edge indices edge out of of k according to the adjacency matrix
    if(length(edges_from_k) != 0) {
        for(a in 1:length(edges_from_k)) {
            for(b in 1:length(edges_from_k)) {
                if(b>a) {
                    e1 <- edges_from_k[a]
                    e2 <- edges_from_k[b]
                    m1 <- edge_models[[e1]]
                    m2 <- edge_models[[e2]]
                    if(is.null(m1) || is.null(m2)) {
                        next
                    }
                    distance <- compare_gmm_models(m1, m2, number_of_gmm_populations)

                    # Just bookkeeping code to record the node names in a standardized order,
                    # together with the name of the 'pivot' node for the adjacent edge pair.
                    four_endpoints <- c(edges[e1, ], edges[e2, ])
                    three_endpoints <- unique(c(edges[e1, ], edges[e2, ]))
                    pivot <- four_endpoints[duplicated(four_endpoints)]
                    endpoints <- setdiff(three_endpoints, pivot)
                    endpoints_s <- sort(c(endpoints[[1]], endpoints[[2]]))

                    stats_k[[count]] <- t(data.frame(c(e1, e2, distance, endpoints_s[1], endpoints_s[2], pivot)))
                    colnames(stats_k[[count]]) = c("edge1", "edge2", "distance", "endpoint1", "endpoint2", "pivot")
                    rownames(stats_k[[count]]) <- c()
                    count <- count + 1
                }
            }
        }
    }
    l <- length(stats_k)
    df_k <- data.frame( edge1=c(rep(1, l)),
                        edge2=c(rep(1, l)),
                        distance=c(rep(0.0, l)),
                        endpoint1=c(rep("X", l)),
                        endpoint2=c(rep("X", l)),
                        pivot=c(rep("X", l)), stringsAsFactors=FALSE)
    colnames(df_k)=c("edge1", "edge2", "distance", "endpoint1", "endpoint2", "pivot")
    if(l != 0) {
        for(i in 1:l) {
            for(j in 1:6) {
                df_k[i, j] <- stats_k[[i]][j]
            }
        }
    }
    df_k$edge1 <- as.numeric(df_k$edge1)
    df_k$edge2 <- as.numeric(df_k$edge2)
    df_k$distance <- as.numeric(df_k$distance)
    cleaned <- na.omit(df_k)
    return(cleaned)
}

#' Virtual edge weight calculation
#' 
#' For each pair of two nodes connected by some triangle, finds the average GMT (W2+GMM) distance among all such triangles.
#' @param eps The edge pair statistics data frame as returned by \code{calculate_edge_pair_statistics}
#' @param log_place File/connection for messages.
#' @return A big 'ol data frame of 'virtual' edges, with \code{source}, \code{target}, and \code{average distance} columns
calculate_average_distance_over_triangles <- function(eps, log_place=stdout()) {
    # *Replace with factor-based version
    eps$average_distance <- c(rep(-1,nrow(eps)))
    eps$distance <- as.numeric(eps$distance)
    tally <- eps[1,"distance"]
    count <- 1
    end1 <- eps$endpoint1
    end2 <- eps$endpoint2
    dists <- eps$distance
    cat("\n", file=log_place, append=TRUE)
    for(i in 2:nrow(eps)) {
        if( (end1[i-1] == end1[i]) && (end2[i-1] == end2[i]) ) {
            tally <- tally + dists[i]
            count <- count + 1
        } else {
            eps[i-1, "average_distance"] <- tally/count
            tally <- dists[i]
            count <- 1
        }
        cat(paste0("\r      ", i, " of ", nrow(eps), " edge pairs collated."), file=log_place, append=TRUE)
    }
    cat("\n", file=log_place, append=TRUE)
    eps[nrow(eps),"average_distance"] <- tally/count

    virtual_edges <- eps[eps$average_distance != -1, c("endpoint1","endpoint2","average_distance")]
    colnames(virtual_edges) <- c("source","target","average_distance")
    virtual_edges <- virtual_edges[order(as.numeric(virtual_edges[, "average_distance"])), ]
    rownames(virtual_edges) <- c()
    return(virtual_edges)
}

make_graph_from_linkage <- function(linkage, node_names, label="linkage"){
    merges <- linkage$merge
    heights <- linkage$height
    last_height <- heights[length(heights)]
    threshold <- length(heights)
    if(last_height == heights[length(heights)-1]) {
        threshold <- which(heights == last_height)[1]-1
    }

    absorption_times <- list()
    k <- 1

    membership_edges <- data.frame(matrix(c(""), nrow=2*threshold, ncol=4), stringsAsFactors=FALSE)
    colnames(membership_edges) <- c("source", "target", label, "type")
    for(i in 1:threshold) {
        membership_edges[2*i-1, 1] <- linkage_index_to_name(merges[i,1], node_names)
        membership_edges[2*i-1, 2] <- paste0("SN",i)
        membership_edges[2*i-1, 3] <- label
        membership_edges[2*i, 1] <- linkage_index_to_name(merges[i,2], node_names)
        membership_edges[2*i, 2] <- paste0("SN",i)
        membership_edges[2*i, 3] <- label

        if(merges[i,1] < 0) {
            absorption_times[[k]] <- c(linkage_index_to_name(merges[i,1], node_names), k)
            k <- k + 1 
        }
        if(merges[i,2] < 0) {
            absorption_times[[k]] <- c(linkage_index_to_name(merges[i,2], node_names), k)
            k <- k + 1 
        }
    }

    new_names <- c(node_names, paste0("SN", c(1:threshold)))
    hierarchy <- graph_from_data_frame(membership_edges, directed=TRUE, vertices = as.data.frame(new_names))
    V(hierarchy)$type <- c(rep("original", length(node_names)), rep("supernode", threshold))
 
    V(hierarchy)$absorption_time <- c(rep(-1*length(absorption_times),length(V(hierarchy))))
    for(j in 1:length(absorption_times)) {
        fname <- absorption_times[[j]][1]
        val <- absorption_times[[j]][2]
        V(hierarchy)[fname]$absorption_time <- -1*as.numeric(val)
    }
 
    V(hierarchy)$int_type <- sapply(V(hierarchy)$type, function(x){if(x=="original"){return(1)} else {return(0)}})
    V(hierarchy)$simple_name <- sapply(V(hierarchy)$name, function(x){if(grepl("SN\\d+",x) ){return("")} else {return(x)}})

    if(length(node_names) != (dim(merges)[1])+1) {
        warning(paste0("      Node number ", length(V(hierarchy)), ", but linkage has ", dim(merges)[1], " merge entries.\n"))        
    }
    return(hierarchy)
}



