# Contact: Jimmy  mathewj2@mskcc.org
#
#' @import mclust
#' @import emdist

suppressMessages(library(mclust)) # Are these calls needed/used?
library(emdist)

#' Compare GMMs with 2-Wasserstein distance
#' 
#' Calculates the W2 distance between Gaussian mixture models, as in Chen et al. 2019.
#' @param m1 First model (as returned by \code{gmm_model_an_edge})
#' @param m2 Second model (as returned by \code{gmm_model_an_edge})
#' @param number_of_gmm_populations Number of populations for the Gaussian Mixture Modeling
#' @return The distance
#' @export
compare_gmm_models <- function(m1, m2, number_of_gmm_populations = 3) {
    probs1 <- m1$probabilities
    probs2 <- m2$probabilities

    if(is.null(probs1) || is.null(probs2)) {
        return(NA)
    }
    if(is.na(probs1) || is.na(probs2)) {
        return(NA)
    }

    n <- number_of_gmm_populations

    # A, B are probability distributions in virtual, small discrete space.
    # The 1:n is the range of indices of the discrete location.
    # These are allowed to be any indices, used later by 'emdr' when looks up cost values to go from one indexed location to another.
    # Non-trivial usage of these indices appears in the Ollivier-Ricci curvature function (possibly included in this package).
    A <- matrix(c(probs1, 1:n), nrow=n, ncol=2)
    B <- matrix(c(probs2, 1:n), nrow=n, ncol=2)

    W2 <- matrix(0, nrow=n, ncol=n)
    # Technically I should be checking that the s1 and s2 matrices below are actually invertible...
    for(x in 1:n) {
        for(y in 1:n) {
            sqrt_plus <- m1$sigma[, ,x]^(1/2)
            W2[x, y] <- norm(as.matrix(m1$mean[, x] - m2$mean[, y]), type="F")^2 +
                        sum(diag(m1$sigma[, ,x] + m2$sigma[, ,y] - 2 * (sqrt_plus * m2$sigma[, ,y] * sqrt_plus)^(1/2) ))
        }
    }
    D <- function(x, y) { sqrt(W2[x, y]) }
    dist <- -1
    tryCatch( {dist <- emdr(A, B, dist=D)},
            error=function(err) {print("emd error")})
    return(dist)
}

#' Gaussian Mixture Modeling of data along an edge
#' 
#' Uses \code{mclust} to model the bivariate joint distribution of data values at endpoint nodes for the given edge.
#' @param edge Vector of endpoint node names
#' @param data The whole node data frame
#' @param number_of_gmm_populations Parameter for \code{mclust}
#' @return A named list of \code{probabilities}, \code{mean}, and \code{sigma}.
#' @export
gmm_model_an_edge <- function(edge, data, number_of_gmm_populations = 3) {
    node_pair_samples <- data[, c(edge[1], edge[2])]
    n <- number_of_gmm_populations
    ps <- NA
    tryCatch(
        {
            mixture<-Mclust(node_pair_samples, G=n, modelNames="VVV")
            if(is.null(mixture)) {
                #cat(paste0("       Mclust model is null for (", edge[1], ", ", edge[2],")\n"))
                return(NULL)
            }
            ps <- mixture$parameters
            return(list(probabilities = ps$pro, mean = ps$mean, sigma = ps$variance$sigma))
        },
        error=function(err){
            cat(paste0("       Modeling failed for (", edge[1], ", ", edge[2], ")"))
            return(NULL)
        })
}

