#' @title 1 test and 1 population binomial mixture model
#' @description 1 test and 1 population binomial mixture model to estimate prevalence and diagnostic test related meassures.
#' @param dataset \code{\link{list}} with the population size and the number of positives. The names of these values must be "pop_size" and "positives" respectively.
#' @param inits \code{\link{list}} with initial conditions for chains. \code{inits} must define initial values for the true prevalence within infected herds, the sensitivity, the specificity and the herd prevalence. The names of these values must be "true_prev_wph", "se", "sp" and "prev_h" respectively.
#' @param priors \code{\link{vector}} with the parameters a and b (Beta distribution) for the true prevalence within infected herds, the sensitivity and the specificity; and with the parameter p (Bernoulli distribution) for the herd prevalence. The names of these values must be: "true_prev_wph_a", "true_prev_wph_b", "se_a", "se_b", "sp_a", "sp_b" and "prev_h".
#' @param pars character vector giving the names of parameters to be monitored. It is passed to the \code{variable.names} argument of the \code{\link{coda.samples}} function. In addition to the parameters specified in \code{inits}, the apparent prevalence, the true prevalence, and the positive and negative predictive values can be monitored if specified in \code{pars} as "app_prev", "true_prev", "ppv" and "npv", respectively.
#' @param n.chains the number of parallel chains for the model. It is passed to the \code{n.chains} argument of the \code{\link{jags.model}} function.
#' @param burn.in the number of iteration to be discarded. It is passed to the \code{n.iter} argument of the \code{\link{update.jags}} function.
#' @param thin thinning interval for monitors. It is passed to the \code{thin} argument of the \code{\link{coda.samples}} function.
#' @param n.iter number of iterations to monitor. It is passed to the \code{n.iter} argument of the \code{\link{coda.samples}} function.
#' #' @details This function creates a text file with the model and it is saved in the working directory.
#' @return A \code{\link{list}} of class \code{mcmc.list}.
#' @references https://dl.dropboxusercontent.com/u/49022/diagnostictests/index.html
#' @export
#' @examples 
#' # Data (initial values for chains automatically generated).
#' dataset <- list(pop_size = 91, positives = 1)
#' 
#' # Priors
#' priors <- list(true_prev_wph_a = 1.8, true_prev_wph_b = 26.74, se_a = 6.28,
#'                se_b = 13.32, sp_a = 212.12, sp_b = 3.13, prev_h = 0.1)
#'                
#' 
#' # Prevalence estimates
#' prev_est <- OneTestOnePopBMM(dataset = dataset, priors = priors, n_iter = 3e3,
#'                              pars = c('true_prev', 'true_prev_wph', 'prev_h'))
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
OneTestOnePopBMM <- function(dataset, inits, priors, pars, n_iter = 1e4, n_chains = 3, burn_in = 1e3, thin = 1) {
    model <- 
        paste(c('model {',
                'positives ~ dbin(app_prev, pop_size)',
                'app_prev <-true_prev * se + (1 - true_prev) * (1 - sp)',
                'ppv <- true_prev * se / app_prev',
                'npv <- (1 - true_prev) * sp / (1 - app_prev)',
                paste0('prev_h ~ dbern(', priors[["prev_h"]], ')'),
                paste0('se ~ dbeta(', priors[["se_a"]], ', ',
                       priors[["se_b"]], ')'),
                paste0('sp ~ dbeta(', priors[["sp_a"]], ', ',
                       priors[["sp_b"]], ')'),
                paste0('true_prev_wph ~ dbeta(',
                       priors[['true_prev_wph_a']],
                       ', ', priors[['true_prev_wph_b']], ')'),
                'true_prev <- prev_h * true_prev_wph',
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