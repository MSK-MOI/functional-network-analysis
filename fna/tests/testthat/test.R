context("Basic functionality testing for network reduction code")
library(grapeplots)

test_that("GMM modeling produces models with meaningful probabilities, means, and variances", {
    cat("\n")
    data <- matrix(0.01*sample(100, 500, replace=TRUE), ncol=5, nrow=100)
    colnames(data) <- paste0("node", c(1:5))
    edge <- c("node2", "node4")

    sink("~/.tmp.saux") # This prevents messages from disrupting the test flow
    gmm <- gmm_model_an_edge(edge, data, number_of_gmm_populations = 3)
    sink(NULL)

    expect_is(gmm$probabilities, "numeric")
    expect_equal(length(gmm$probabilities), 3)

    expect_is(gmm$mean, "matrix")
    expect_equal(dim(gmm$mean), c(2,3))
    expect_true(is.numeric(gmm$mean))

    expect_is(gmm$sigma[,,1], "matrix")
    expect_equal(dim(gmm$sigma[,,1]), c(2,2))
    expect_true(is.numeric(gmm$sigma[,,1]))
})
