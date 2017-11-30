#' @title 1 test and 1 population binomial model
#' @description 1 test and 1 population binomial model to estimate prevalence and diagnostic test related meassures.
#' @param dataset \code{\link{list}} with the population size and the number of positives. The names of these values must be "pop_size" and "positives" respectively.
#' @param inits \code{\link{list}} with initial conditions for chains. \code{inits} must define initial values for the true prevalence, the sensitivity and the specificity. The names of these values must be "true_prev", "se" and "sp" respectively.
#' @param priors \code{\link{vector}} with the parameters a and b (Beta distribution) for the true prevalence, the sensitivity and the specificity. The names of these values must be: "true_prev_a", "true_prev_b", "se_a", "se_b", "sp_a" and "sp_b".
#' @param pars character vector giving the names of parameters to be monitored. It is passed to the \code{variable.names} argument of the \code{\link{coda.samples}} function. In addition to the parameters specified in \code{inits}, the apparent prevalence and the positive and negative predictive values can be monitored if specified in \code{pars} as "app_prev", "ppv" and "npv", respectively.
#' @param n_chains the number of parallel chains for the model. It is passed to the \code{n.chains} argument of the \code{\link{jags.model}} function.
#' @param burn_in the number of iteration to be discarded. It is passed to the \code{n.iter} argument of the \code{\link{update.jags}} function.
#' @param thin thinning interval for monitors. It is passed to the \code{thin} argument of the \code{\link{coda.samples}} function.
#' @param n_iter number of iterations to monitor. It is passed to the \code{n.iter} argument of the \code{\link{coda.samples}} function.
#' @details This function creates a text file with the model and it is saved in the working directory.
#' @return A \code{\link{list}} of class \code{mcmc.list}.
#' @references https://dl.dropboxusercontent.com/u/49022/diagnostictests/index.html
#' @export
#' @examples 
#' # Data.
#' dataset <- list(pop_size = 91, positives = 1)
#' 
#' # Initial conditions for chains.
#' inits <- list(list(true_prev = 0.05, se = 0.8, sp = 0.9),
#'               list(true_prev = 0.02, se = 0.3, sp = 0.7),
#'               list(true_prev = 0.09, se = 0.1, sp = 0.5))
#' 
#' # Priors.
#' priors <- c(true_prev_a = 1, true_prev_b = 1,
#'             se_a = 6.28, se_b = 13.32, sp_a = 212.12, sp_b = 3.13)
#' 
#' # Prevalence estimate.
#' prev_est <- OneTestOnePopBM(dataset = dataset, inits = inits, n_iter = 3e3,
#'                             priors = priors, pars = 'true_prev')
#' 
#' summary(prev_est)
#' 
#' # Diagnostic plots.
#' library(coda); library(ggmcmc)
#' gelman.diag(prev_est)
#' gelman.plot(prev_est)
#' gg_res <- ggs(prev_est)
#' ggs_traceplot(gg_res)
#' ggs_density(gg_res)
#' ggs_histogram(gg_res, bins = 100)
#' ggs_compare_partial(gg_res)
#' ggs_running(gg_res)
#' ggs_autocorrelation(gg_res)
OneTestOnePopBM <- function(dataset, inits, priors, pars, n_iter = 1e4,  n_chains = 3, burn_in = 1e3, thin = 1) {
    model <- 
        paste(c('model {',
                'positives ~ dbin(app_prev, pop_size)',
                'app_prev <- true_prev * se + (1 - true_prev) * (1 - sp)',
                'ppv <- true_prev * se / app_prev',
                'npv <- (1 - true_prev) * sp / (1 - app_prev)',
                paste0('se ~ dbeta(', priors[["se_a"]], ', ',
                       priors[["se_b"]], ')'),
                paste0('sp ~ dbeta(', priors[["sp_a"]], ', ',
                       priors[["sp_b"]], ')'),
                paste0('true_prev ~ dbeta(', priors[["true_prev_a"]], ', ',
                       priors[["true_prev_b"]], ')'),
                '}'),
              collapse = '\n')  
    writeLines(model, 'model.txt')
    jags_mod <- jags.model(file = 'model.txt',
                           data = dataset,
                           inits = inits, 
                           n.chains = n_chains)
    update(jags_mod, n.iter = burn_in)
    coda.samples(jags_mod,
                 variable.names = pars,
                 n.iter = n_iter,
                 thin = thin)
}