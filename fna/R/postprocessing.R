# Contact: Jimmy  mathewj2@mskcc.org
#
#' @import stringr
#' @import igraph
#' @import rols
#' @importFrom stats median
#' @importFrom utils head
#' @importFrom utils write.table


library(stringr)
suppressMessages(library(igraph))
suppressMessages(library(rols)) # Gotten from BiocManager::install("rols")

#' Basic annotation of network reduction
#' 
#' Annotate a network reduction run with feature clusters and with profiles of data sets along sample clusters. Feature names of supplied data set should at least somewhat match the node names of the supplied graph.
#' @param graph igraph as returned by \link{generate_reduction}.
#' @param node_data_file Data file with numerical node sample data, with node names. Can be simple CSV, or GCT format. Does not need to be the original data used in calculating a network reduction!
#' @param cluster_file 2-column CSV file describing a sample clustering.
#' @param correlation_transformation (optional) (see \link{generate_reduction}).
#' @param normalization Normalization for the data set (see \link{generate_reduction}).
#' @param verbose (optional) Whether to print messages to console. Default FALSE.
#' @param log_file (optional) A file to log messages. Default NA.
#' @return \code{graph}
#' @export
annotate_with_clusters <- function(graph,
                                   node_data_file,
                                   cluster_file,
                                   correlation_transformation=c("none", "pearson", "spearman"),
                                   normalization=c("dividebymean", "zscore", "none"),
                                   verbose=FALSE,
                                   log_file=NA
                                   ) {
    correlation_transformation <- match.arg(correlation_transformation)
    normalization <- match.arg(normalization)

    data <- get_sample_data(node_data_file)
    if(correlation_transformation != "none") {
        data <- cor(data, method=correlation_transformation)
    } else {
        data <- normalize(data, normalization_method=normalization)
    }

    if(!is.na(cluster_file)) {
        descriptor <- sub("\\.(csv|CSV)$", "", basename(cluster_file))
        assignments <- read.csv(cluster_file, stringsAsFactors=FALSE)
        rownames(assignments) <- assignments[,1]
        colnames(assignments) <- c("id", "cluster")

        all_ids <- intersect(assignments[,"id"], rownames(data))
        all_features <- intersect(V(graph)$name, colnames(data))
        data_local <- normalize(data[all_ids, all_features], normalization_method="zscore")
        if(any(is.na(data_local))) {
            warning("Some NA values in the data set along the nodes of the graph.")
            return(graph)
        }

        profiles <- list()
        fac <- as.factor(assignments[,"cluster"])
        log_message(paste0("Annotating based on ", length(all_ids), " samples (of ", dim(data)[1]," total possible) in groups of sizes:\n", paste(as.vector(table(fac)), collapse=" "), "\n"), verbose=verbose, log_file=log_file)
        for(i in 1:length(fac)) {
            level_i <- levels(fac)[i]
            ids_i <- intersect(assignments[assignments[,"cluster"]==level_i, "id"], all_ids)
            if(any(is.na(data_local[ids_i,]))) {
                warning(paste0("Some NA values in the data of level ", level_i, "."))
                return(graph)
            }
            profile_in_data <- colMeans(data_local[ids_i,])
            if(any(is.na(profile_in_data))) {
                warning(paste0("Some NA values in the profiles of level ", level_i, "."))
                return(graph)
            }
            profile_on_nodes <- sapply(V(graph)$name, FUN=function(x) {return(profile_in_data[x])})
            # profile_on_nodes is allowed to have NA values, since not all nodes are necessarily in the data set.
            vertex_attr(graph, paste0(descriptor, "_", level_i)) <- profile_on_nodes
        }
    }
    return(graph)
}

#' Annotate network reduction with feature clusters
#' 
#' Annotate a network reduction run with labels for groups of nodes/features.
#' @param graphml_file File name for GraphML formatted graph.
#' @param feature_cluster_file 2-column CSV file describing a feature clustering.
#' @param outdir Subdirectory for output (default NA, meaning current directory).
#' @export
annotate_graphml_with_feature_clusters <- function(graphml_file, feature_cluster_file, outdir=NA) {
    assignments <- read.csv(feature_cluster_file, stringsAsFactors=FALSE)
    rownames(assignments) <- assignments[,1]
    descriptor <- colnames(assignments)[2]
    
    graph <- read_graph(graphml_file, format="graphml")

    all_assignments <- sapply(V(graph)$name,
        function(x){
            if(x %in% rownames(assignments)) {
                return(assignments[x, 2])
            } else {
                return(-1)
            }
        })
    vertex_attr(graph, descriptor) <- all_assignments

    outfile <- paste0(sub("\\.(graphml|GRAPHML)$", "", graphml_file), "_annotated.graphml")
    if(!is.na(outdir)) {
        outfile <- paste(outdir, basename(outfile), sep="/")
    }
    write_graph(graph, outfile, format="graphml")
}

#' Annotate a GraphML file using Gene Annotation File database
#' 
#' Annotates a GraphML file whose nodes are genes. See \link{gaf_annotate_igraph}.
#' @param graphml_file GraphML file with 'name' attribute with valid HUGO gene symbols.
#' @param gaf_file GAF (Gene Annotation Format) file (default \code{goa_human.gaf}) to get annotations from.
#' @param outdir Subdirectory for output (default NA, meaning current directory).
#' @param rename If FALSE, does not give a new name to the output file (possibly overwriting the input GraphML file). If TRUE, renames with "_annotated" appended to filename (default TRUE).
#' @param linelimit (optional) Limit on number of lines of gaf_file to read. Default NA (no limit).
#' @param verbose (optional) Whether to print messages to console. Default FALSE.
#' @param log_file (optional) A file to log messages. Default NA.
#' @export
gaf_annotate_graphml <- function(graphml_file, gaf_file="goa_human.gaf", outdir=NA, rename=TRUE, linelimit=NA, verbose=FALSE, log_file=NA) {
    #warning("igraph library's GraphML parser doesn't always load edge attributes from GraphML.")
    if(!grepl("\\.(graphml|GRAPHML)$", graphml_file)) {
        warning("graphml_file not recognized as being in GraphML file format.")
        return()
    }

    graph <- read_graph(graphml_file, format = "graphml")
    # graph <- igraph::simplify(read_graph(graphml_file, format = "graphml"))
    graph <- gaf_annotate_igraph(graph, gaf_file=gaf_file, linelimit=linelimit)   #egg

    # appropriate time for this check?
    if(length(intersect(c("name", "absorption_time"), list.vertex.attributes(graph))) != 2) {
        warning("Graph doesn't have name or doesn't have absorption_time attributes.")
        # return()
    }

    if(!rename) {
        outfile <- graphml_file
    } else {
        outfile <- paste0(sub("\\.(graphml|GRAPHML)$", "", graphml_file), "_annotated.graphml")
    }

    if(!is.na(outdir)) {
        outfile <- paste(outdir, basename(outfile), sep="/")
    }

    write_graph(graph, outfile, format="graphml")
}

#' Annotate an igraph using Gene Annotation File database
#' 
#' Annotates an igraph object whose nodes are labelled with a 'name' field with HUGO gene symbols, using the annotations in a GAF file. Defaults to a file downloaded from the EBI (\url{https://www.ebi.ac.uk/GOA/downloads}). To avoid the slowdown of repeated and redundant REST API calls, a local cache is maintained. For reuse, it is written to file 'go_terms_local.cache'.
#' @param graph igraph with a 'name' node attribute with valid HUGO gene symbols.
#' @param gaf_file GAF (Gene Annotation Format) file (default \code{goa_human.gaf}) to get annotations from.
#' @param offline (optional) If FALSE, will use package 'rols' to make web queries to lookup GO terms. If TRUE will only use local cache file 'go_terms_local.cache' (default FALSE).
#' @param linelimit (optional) Limit on number of lines of gaf_file to read. Default NA (no limit).
#' @param verbose (optional) Whether to print messages to console. Default FALSE.
#' @param log_file (optional) A file to log messages. Default NA.
#' @export
gaf_annotate_igraph <- function(graph, gaf_file="goa_human.gaf", offline=FALSE, linelimit=NA, verbose=FALSE, log_file=NA) {
    if(!grepl("\\.(gaf|GAF)$", gaf_file)) {
        warning("gaf_file not recognized as being in GAF file format.")
        return()
    }

    if(length(intersect("name", list.vertex.attributes(graph))) != 1) {
        warning("igraph object does not have a 'name' node attribute. Trying 'label'.")
        if(length(intersect("label", list.vertex.attributes(graph))) != 1) {
            warning("igraph object doesn't have 'label' node attribute either.")
            return()
        }
    }

    pattern <- get_gaf_line_pattern()
    go <- Ontology("go")
    if(file.exists("go_terms_local.cache")) {
        lookuptable <- readRDS("go_terms_local.cache")
        lookedup_label <- lookuptable[[1]]
        lookedup_description <- lookuptable[[2]]
    } else {
        lookedup_label <- c()
        lookedup_description <- c()
        if(offline) {
            stop("Selected offline mode for annotation, but there is no local cache file.")
        }
    }

    # This function does the reading/parsing of a GAF file line, doesn't refer to the input graph, but I think there was some reason I kept it local to the 'gaf_annotate_igraph' function. Some scope issue.
    process_line <- function(line, pattern) {
        matches <- str_match(line, pattern)

        if(length(matches) < 16 || any(is.na(matches))) {
            print(matches)
            warning(paste0("Tried to process a bad line: ", line))
            return(NA)
        }

        hugo_name <- matches[4]
        go_entity_for_involvement <- matches[6]
        pubmed_id <- matches[7]
        evidence_code <- matches[8]
        aspect <- matches[10]
        full_name <- matches[11]

        all_aspects <- data.frame(gocat=c("process","function","component"), row.names=c("P","F","C"), stringsAsFactors=FALSE)
        aspect <- all_aspects[aspect,]

        lookup <- ols_call(go_entity_for_involvement)
        if(length(lookup)==1 && lookup=="offline") {
            return("offline")
        }

        annotation <- lookup[[1]][1]
        definition <- lookup[[1]][2]
        lookedup_label <- lookup[[2]]
        lookedup_description <- lookup[[3]]
        return(list(c(hugo_name, full_name, annotation, definition, aspect, pubmed_id, evidence_code), lookedup_label, lookedup_description))
    }

    # Same comment applies as for the 'process_line' function above
    ols_call <- function(go_entity_for_involvement) {
        if(go_entity_for_involvement %in% names(lookedup_label)) {
            return(list(c(lookedup_label[go_entity_for_involvement], lookedup_description[go_entity_for_involvement]), lookedup_label, lookedup_description))
        } else {
            if(offline) {
                return("offline")
            }
            t <- term(go, go_entity_for_involvement)
            annotation <- termLabel(t)
            definition <- termDesc(t)
            # Note that 'lookedup_...'' are in the scope here since we're doing fancy nested functions. Since R doesn't allow passing objects by reference for convenient in-place modification. Same applies to 'go'.
            if(is.null(annotation)) {
                warning(paste0("Annotation for ", go_entity_for_involvement, " was null."))
                annotation <- "ANNOTATION WAS NULL"
            } else if(is.na(annotation)) {
                warning(paste0("Annotation for ", go_entity_for_involvement, " was not found."))
                annotation <- "NO ANNOTATION FOUND"
            }

            if(is.null(definition)) {
                warning(paste0("Definition for ", go_entity_for_involvement, " was null."))
                definition <- "DEFINITION WAS NULL"
            }
            if(is.na(definition)) {
                warning(paste0("Definition for ", go_entity_for_involvement, " was not found."))
                definition <- "NO DEFINITION FOUND"
            }

            lookedup_label <- append(lookedup_label, annotation)
            lookedup_description <- append(lookedup_description, definition)
            names(lookedup_label)[length(lookedup_label)] <- go_entity_for_involvement
            names(lookedup_description)[length(lookedup_description)] <- go_entity_for_involvement
            return(list(c(annotation, definition), lookedup_label, lookedup_description))
        }
    }

    number_lines <- as.integer(str_extract(system(paste("wc", "-l", gaf_file), intern=TRUE), "\\d+"))
    log_message(paste0(number_lines, " lines to scan."), verbose=verbose, log_file=log_file)

    l <- c(rep(NA, number_lines))
    hugo_name <- l
    full_name <- l
    annotation <- l
    definition <- l
    aspect <- l
    pubmed_id <- l
    evidence_code <- l

    if(length(intersect(c("name"), list.vertex.attributes(graph))) != 1) {
        # Trying 'label'
        if(length(intersect(c("label"), list.vertex.attributes(graph))) == 1) {
            vertex_attr(graph, "name", V(graph)) <- V(graph)$label
        } else {
            stop("No 'label' node attribute.")
        }
    }

    names <- V(graph)$name
    connection <- file(gaf_file, open="r")
    k <- 1
    i <- 1

    log_message("\n", verbose=verbose, log_file=log_file)
    start_time <- Sys.time()
    while(TRUE) {
        if(!is.na(linelimit)) {
            if(i == linelimit) {
                break
            }
        }

        i <- i + 1

        stop_time <- Sys.time()
        elapsed <- format(round(difftime(stop_time, start_time, 2, units="mins"), 2), nsmall = 2)
        estimate <- format(round((number_lines/i)*difftime(stop_time, start_time, units="mins"), 2), nsmall=2) 
        if(verbose) {
            cat(paste0(k-1, " relevant lines found in ", gaf_file, " out of ", i, " total lines read. (Elapsed ", elapsed, " of expected total ", estimate, ")\r"))
        }

        line <- readLines(connection, n=1, warn=FALSE)
        if(is.na(line) || length(line)==0) {
            break
        }
        first_char <- substr(line, 1, 1)
        if(first_char == "!") {
            next
        }

        fields <- strsplit(line, "\\t")[[1]]

        if(length(fields) < 15) {
            warning(paste0("Not enough complete fields in line. Assuming end of data."))
            break
        }

        if(!(fields[3] %in% names)) {
            next
        }
        
        record <- process_line(line, pattern)
        if(length(record)==1 && record=="offline") {
            next
        }
        for(w in 1:length(record)) {
            if(any(is.na(record[[w]]))) {
                warning(paste0("NA value in position ", w, " line ", i, " and item ", k))
            }
        }
        lookedup_label <- record[[2]]
        lookedup_description <- record[[3]]
        hugo_name[k] <- record[[1]][1]
        full_name[k] <- record[[1]][2]
        annotation[k] <- record[[1]][3]
        definition[k] <- record[[1]][4]
        aspect[k] <- record[[1]][5]
        pubmed_id[k] <- record[[1]][6]
        evidence_code[k] <- record[[1]][7]
        k <- k + 1
    }
    close(connection)
    log_message("\n", verbose=verbose, log_file=log_file)

    if(!offline) {
        saveRDS(list(lookedup_label, lookedup_description), "go_terms_local.cache")
    }
    df <- data.frame(hugo_name=hugo_name, full_name=full_name, annotation=annotation, definition=definition, aspect=aspect, pubmed_id=pubmed_id, evidence_code=evidence_code, stringsAsFactors=FALSE)
    df <- df[1:(k-1),]

    full_names <- df$full_name[!duplicated(df$hugo_name)]
    names(full_names) <- df$hugo_name[!duplicated(df$hugo_name)]
    V(graph)$full_name <- full_names[V(graph)$name]

    dfr <- list(C=df[df$aspect == "component",], F=df[df$aspect == "function",], P=df[df$aspect == "process",])
    aspects <- c("C","F","P")
    tb <- list()

    for(k in 1:3) {
        dfk <- dfr[[k]]
        max_levels <- 100

        # Ranking mechanism. Calculate average between-node distance on hierarchy graph for subset with given annotation.
        leaf_nodes <- V(graph)$name
        # Definitely only works on the hierarchy version now
        edge_attr(graph, "weight", E(graph)) <- -1*(V(graph)[ends(graph,E(graph))[,2]])$absorption_time_refined
        node_distances <- distances(graph, v=leaf_nodes, to=leaf_nodes, mode="all")

        annotations <- unique(dfk$annotation)
        medians_and_pvalues <- sapply(annotations,
            function(annotation){
                annotated_genes <- unique(dfk$hugo_name[dfk$annotation == annotation])
                in_bc <- intersect(rownames(node_distances), annotated_genes)
                if(length(in_bc)<=3) {
                    return(c(Inf, 1.0))
                } else {
                    # return((1.0/length(annotated_genes)) * fancymedian(node_distances[in_bc, in_bc]))
                    # return( fancymedian(node_distances[in_bc, in_bc], trials=1000, nodes=in_bc, node_distances=node_distances))
                    return( fancymedian(node_distances[in_bc, in_bc], trials=1000, nodes=in_bc, node_distances=node_distances, type="average"))
                }
            })
        per_node_cluster_means <- medians_and_pvalues[1,]
        pvalues <- medians_and_pvalues[2,]
        number_genes <- sapply(annotations,
            function(annotation){
                annotated_genes <- unique(dfk$hugo_name[dfk$annotation == annotation])
                return(length(annotated_genes))
            })
        ranks <- rank(per_node_cluster_means, ties.method="first")
        names(ranks) <- annotations
        names(per_node_cluster_means) <- annotations
        names(pvalues) <- annotations
        dfk$ranks <- sapply(dfk$annotation, FUN=function(annotation){return(ranks[annotation])})

        tb[[k]] <- data.frame(average_stat=per_node_cluster_means[1:max_levels], pvalues=pvalues[1:max_levels], number_genes=number_genes[1:max_levels], stringsAsFactors=FALSE)
        rownames(tb[[k]]) <- names(per_node_cluster_means[1:max_levels])
        tb[[k]] <- tb[[k]][order(tb[[k]][,"pvalues"]),]
        log_message(aspects[k], verbose=verbose, log_file=log_file)
        print(head(tb[[k]], n=10))
        write.table(tb[[k]], paste0("log", aspects[k], ".txt"), quote=FALSE, row.names=TRUE)

        # Record annotation assignments in order of rank
        for(i in 1:max_levels) {
            # namesi <- dfk$hugo_name[dfk$ranks == i]
            # namesi <- namesi[!duplicated(namesi)]
            # annotationi <- dfk$annotation[which.max(dfk$ranks == i)]
            annotationi <- rownames(tb[[k]])[i]
            namesi <- dfk$hugo_name[dfk$annotation == annotationi]
            namesfilter <- sapply(V(graph)$name, FUN=function(x){return(x %in% namesi)})

            attribute_index <- (100*(k-1) + (i-1)) + 1
            vertex_attr(graph, paste0("GOA", attribute_index, "::p=", signif(pvalues[annotationi], digits=3), "::", annotationi), V(graph)[namesfilter]) <- "TRUE"
            vertex_attr(graph, paste0("GOA", attribute_index, "::p=", signif(pvalues[annotationi], digits=3), "::", annotationi), V(graph)[!namesfilter]) <- "FALSE"
        }
        log_message("\n\n", verbose=verbose, log_file=log_file)

    }

    return(graph)
}

fancymedian <- function(mm, trials=NA, nodes=NA, node_distances=NA, type=c("average")) {
    vv <- as.vector(mm)
    vv <- vv[(!is.na(vv)) & (!is.nan(vv)) & !(vv==Inf)]
    vv <- vv[vv > 0]
    if(length(vv) < 6) {
        if(is.na(trials)) {
            return(Inf)
        } else {
            return(c(Inf, 1.0))
        }
    }
    if(type == "median") {
        fm <- median(vv)
    }
    if(type == "average") {
        fm <- mean(vv)
    }
    if(fm == Inf) {
        return(c(Inf, 1.0))
    }
    if(is.na(trials)) {
        return(fm)
    } else {
        number_nodes <- dim(node_distances)[1]
        high_count <- 0
        low_count <- 0
        for(i in c(1:trials)) {
            new_nodes <- sample(c(1:number_nodes), length(nodes), replace=FALSE)
            new_fm <- fancymedian(node_distances[new_nodes, new_nodes])
            if(new_fm == Inf) {
                next
            }
            if(new_fm <= fm) {
                low_count <- low_count + 1
            } else {
                high_count <- high_count + 1
            }
        }
        if(low_count + high_count == 0) {
            return(c(fm, 1.0))
        } else {
            return(c(fm, low_count / (low_count + high_count)))
        }
    }
}


#' Gene Annotation File line pattern
#' 
#' Regex pattern for matching fields in GAF files. Reverse-engineered from databases downloaded from the EBI (goa_human.gaf) in June 2019.
#' @return The regex pattern.
#' @export
get_gaf_line_pattern <- function() {
    #                                           number          example
    #
    # 1    DB                       required    1               UniProtKB
    # 2    DB Object ID             required    1               P12345
    # 3    DB Object Symbol         required    1               PHO3                 protein name
    # 4    Qualifier                optional    0 or greater    NOT
    # 5    GO ID                    required    1               GO:0003993           GO entity the protein is associated with. The key thing for us.
    # 6    (|DB:Reference)          required    1 or greater    PMID:2676709
    # 7    Evidence Code            required    1               IMP
    # 8    With (or) From           optional    0 or greater    GO:0000346
    # 9    Aspect                   required    1               F
    # 10   DB Object Name           optional    0 or 1          Toll-like receptor 4
    # 11   DB Object (|Synonym)     optional    0 or greater    hToll|Tollbooth
    # 12   DB Object Type           required    1               protein
    # 13   Taxon(|taxon)            required    1 or 2          taxon:9606
    # 14   Date                     required    1               20090118
    # 15   Assigned By              required    1               SGD
    # 16   Annotation Extension     !optional!  0 or greater    part_of(CL:0000576)
    # 17   Gene Product Form ID     !optional!  0 or 1          UniProtKB:P12345-2
    p <- c()
    p[1] = ""
    p[2] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<\\: ]+)\\t"
    p[3] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<\\: ]+)\\t" 
    p[4] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<\\:/ ]+)\\t"
    p[5] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<\\:\\|/ ]*)\\t"
    p[6] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<\\:/ ]+)\\t"
    p[7] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\<:/ ]+)\\t"
    p[8] = "(\\w{2,3})\\t"
    p[9] = "([\\w\\d_\\-\\+\\:\\[\\]\\.\\,\\>\\<\\| ]*)\\t"
    p[10] = "([PFC])\\t"
    p[11] = "([\\w\\d_\\-\\+\\,\\(\\)\\'\\,\\>\\<\\./\\[\\] \\:]+)\\t"
    p[12] = "([\\w\\d_\\-\\+\\|/\\[\\]\\.\\,\\>\\<\\'\\* ]*)\\t"
    p[13] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\< /]+)\\t"
    p[14] = "([\\w\\d_\\-\\+\\[\\]:\\|\\.\\,\\>\\< ]*)\\t"
    p[15] = "(\\d{8})\\t"
    p[16] = "([\\w\\d_\\-\\+\\[\\]\\.\\,\\>\\< ]+)\\t"
    p[17] = "([\\w\\d()_\\-\\+\\[\\]\\:\\|\\,\\>\\<\\. ]*)\\t?"
    p[18] = "([\\w\\d_\\-\\+\\[\\]\\:\\.\\,\\>\\< ]*)"

    return(paste0(p, collapse=""))
}


