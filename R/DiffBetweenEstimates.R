#' @title Difference between estimated parameters
#' @description Vectorized difference between chains of two objects of \code{\link{class}} \code{mcmc.list}.
#' @param estimates \code{\link{list}} with two objects of class \code{mcmc.list}. These objects must have the same number of chains and the chains must have the same length.
#' @return A list of class \code{mcmc.list}.
#' @export
#' @examples
#' 
#' # Priors.
#' priors <- c(true_prev_a = 1, true_prev_b = 1,
#'             se_a = 6.28, se_b = 13.32, sp_a = 212.12, sp_b = 3.13)
#'
#' # First estimate.
#' dataset1 <- list(pop_size = 100, positives = 5)
#' prev_est1 <- OneTestOnePopBM(dataset = as.list(dataset1), n_iter = 3e3,
#'                                   priors = priors, pars = "true_prev",
#'                                   burn_in = 5e2)
#'
#' # Second estimate.
#' dataset2 <- list(pop_size = 91, positives = 1)
#' prev_est2 <- OneTestOnePopBM(dataset = as.list(dataset2), n_iter = 3e3,
#'                                    priors = priors, pars = "true_prev",
#'                                    burn_in = 5e2)
#' 
#' # Estimated difference.
#' diffs <- DiffBetweenEstimates(list(prev_est1, prev_est2))
#' summary(diffs)
#' 
#' # Diagnostic plots.
#' library(coda); library(ggmcmc)
#' gelman.diag(diffs)
#' gelman.plot(diffs)
#' gg_res <- ggs(diffs)
#' ggs_traceplot(gg_res)
#' ggs_density(gg_res)
#' ggs_histogram(gg_res, bins = 100)
#' ggs_compare_partial(gg_res)
#' ggs_running(gg_res)
#' ggs_autocorrelation(gg_res)
DiffBetweenEstimates <- function(estimates) {
    diffs <- estimates[[1]]
#     temp <- attr(diffs[[1]], "dimnames")[[2]]
#     temp <- paste0(temp, "1 - ", temp, "2  ")
#     for (i in 1:length(diffs)) {
#         attr(diffs[[i]], "dimnames")[[2]] <- temp
#     }
    atts <- attributes(diffs[[1]])
    atts$dim <- c(dim(diffs[[1]])[1], dim(diffs[[1]])[2] * 2)
    par_names <- atts$dimnames[[2]]
    temp.names1 <- temp.names2 <- c()
    for (i in 1:length(par_names)) {
        temp.names1[i] <- paste0(par_names[i], "1 - ", par_names[i], "2")
        temp.names2[i] <- paste0("Pr(", par_names[i], "1 > ", par_names[i], "2)")
    }
    atts$dimnames[[2]] <- c(temp.names1, temp.names2)
    ProportionOfPositives <- function(x) {
        # step function in jags.
        x[x > 0] <- 1
        x[x <= 0] <- 0
        x}
    for (i in 1:length(diffs)) {
        diffs[[i]] <- diffs[[i]] - estimates[[2]][[i]]
        temp <- apply(diffs[[i]], 2, ProportionOfPositives)
        diffs[[i]] <- cbind(diffs[[i]], temp)
        attributes(diffs[[i]]) <- atts
    }
    return(diffs)
}
