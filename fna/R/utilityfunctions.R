# Contact: Jimmy  mathewj2@mskcc.org
#
#' @importFrom stats var
library(CePa)
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
        graph <- igraph::simplify(read_graph(edge_list_file, format = "graphml"))
    } else if(grepl("\\.(ncol|NCOL)$", edge_list_file)) {
        graph <- igraph::simplify(read_graph(edge_list_file, format = "ncol", directed = FALSE))
    } else {
        stop(paste0("Graph file (", edge_list_file, ") format not supported (only graphml and ncol supported)."))
    }
    return(graph)
}

get_common_names <- function(data, graph, verbose=FALSE, log_file=NA) {
    node_names <-V(graph)$name
    sample_attribute_names <- colnames(data)

    d1 <- setdiff(node_names, sample_attribute_names)
    d2 <- setdiff(sample_attribute_names, node_names)

    nodes <- intersect(node_names, sample_attribute_names)
    if(length(nodes) <= 3) {
        stop("3 or fewer nodes in the network with data.")
    }
    log_message(paste0("      Number of nodes in network not in the data set: ", length(d1), "\n"), verbose=verbose, log_file=log_file)
    if(length(d1) > 0) {
        log_message(paste0("      ", paste0(d1[1:min(length(d1), 5)], collapse=" ") ), verbose=verbose, log_file=log_file)
        if(length(d1) > 5) {
            log_message(" ...", verbose=verbose, log_file=log_file)
        }
        log_message("\n", verbose=verbose, log_file=log_file)
    }
    log_message(paste0("      Number of nodes in the data set not in the network: ", length(d2), "\n"), verbose=verbose, log_file=log_file)
    if(length(d2) > 0) {
        log_message(paste0("      ", paste0(d2[1:min(length(d2), 5)], collapse=" ") ), verbose=verbose, log_file=log_file)
        if(length(d2) > 5) {
            log_message(" ...", verbose=verbose, log_file=log_file)
        }
        log_message("\n", verbose=verbose, log_file=log_file)
    }

    log_message(paste0("      Number of nodes in common: ", length(nodes), "\n"), verbose=verbose, log_file=log_file)
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

log_message <- function(message, verbose=FALSE, log_file = NA) {
    if(verbose) {
        cat(message)
    } else {
        if(!is.na(log_file)) {
            write(message, file=log_file, append=TRUE)
        }
    }
}
