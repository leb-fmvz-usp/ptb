#' @title Difference between estimated parameters
#' @description Vectorized difference between chains of two objects of \code{\link{class}} \code{mcmc.list}.
#' @param estimates \code{\link{list}} with two objects of class \code{mcmc.list}. These objects must have the same number of chains and the chains must have the same length.
#' @return A list of class \code{mcmc.list}.
#' @export
#' @examples
#' 
#' # Priors.
#' priors <- c(true.prev.a = 1, true.prev.b = 1,
#'             se.a = 6.28, se.b = 13.32, sp.a = 212.12, sp.b = 3.13)
#'
#' # First estimate.
#' dataset1 <- list(pop.size = 100, positives = 5)
#' prev.est1 <- OneTestOnePopBM(dataset = as.list(dataset1), n.iter = 3e3,
#'                                   priors = priors, pars = 'true.prev',
#'                                   burn.in = 5e2)
#'
#' # Second estimate.
#' dataset2 <- list(pop.size = 91, positives = 1)
#' prev.est2 <- OneTestOnePopBM(dataset = as.list(dataset2), n.iter = 3e3,
#'                                    priors = priors, pars = 'true.prev',
#'                                    burn.in = 5e2)
#' 
#' # Estimated difference.
#' diffs <- DiffBetweenEstimates(list(prev.est1, prev.est2))
#' summary(diffs)
#' 
#' # Diagnostic plots.
#' library(coda); library(ggmcmc)
#' gelman.diag(diffs)
#' gelman.plot(diffs)
#' gg.res <- ggs(diffs)
#' ggs_traceplot(gg.res)
#' ggs_density(gg.res)
#' ggs_histogram(gg.res, bins = 100)
#' ggs_compare_partial(gg.res)
#' ggs_running(gg.res)
#' ggs_autocorrelation(gg.res)
DiffBetweenEstimates <- function(estimates) {
    diffs <- estimates[[1]]
#     temp <- attr(diffs[[1]], "dimnames")[[2]]
#     temp <- paste0(temp, '1 - ', temp, '2  ')
#     for (i in 1:length(diffs)) {
#         attr(diffs[[i]], "dimnames")[[2]] <- temp
#     }
    atts <- attributes(diffs[[1]])
    atts$dim <- c(dim(diffs[[1]])[1], dim(diffs[[1]])[2] * 2)
    par.names <- atts$dimnames[[2]]
    temp.names1 <- temp.names2 <- c()
    for (i in 1:length(par.names)) {
        temp.names1[i] <- paste0(par.names[i], '1 - ', par.names[i], '2')
        temp.names2[i] <- paste0('Pr(', par.names[i], '1 > ', par.names[i], '2)')
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
