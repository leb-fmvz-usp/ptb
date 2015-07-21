#' @title 1 test and 1 population binomial model
#' @description 1 test and 1 population binomial model to estimate prevalence and diagnostic-test-related meassures.
#' @param dataset \code{\link{list}} with the population size and the number of positives. The names of these values must be "pop.size" and "positives" respectively.
#' @param inits \code{\link{list}} with initial conditions for chains. \code{inits} must define initial values for the true prevalence, the sensitivity and the specificity. The names of these values must be "true.prev", "se" and "sp" respectively.
#' @param priors \code{\link{vector}} with the parameters a and b (Beta distribution) for the true prevalence, the sensitivity and the specificity. The names of these values must be: "true.prev.a", "true.prev.b", "se.a", "se.b", "sp.a" and "sp.b".
#' @param pars character vector giving the names of parameters to be monitored. It is passed to the \code{variable.names} argument of the \code{\link{coda.samples}} function. In addition to the parameters specified in \code{inits}, the apparent prevalence and the positive and negative predictive values can be monitored if specified in \code{pars} as "app.prev", "ppv" and "npv", respectively.
#' @param n.chains the number of parallel chains for the model. It is passed to the \code{n.chains} argument of the \code{\link{jags.model}} function.
#' @param burn.in the number of iteration to be discarded. It is passed to the \code{n.iter} argument of the \code{\link{update.jags}} function.
#' @param thin thinning interval for monitors. It is passed to the \code{thin} argument of the \code{\link{coda.samples}} function.
#' @param n.iter number of iterations to monitor. It is passed to the \code{n.iter} argument of the \code{\link{coda.samples}} function.
#' @details This function creates a text file with the model and it is saved in the working directory.
#' @return A \code{\link{list}} of class \code{mcmc.list}.
#' @references https://dl.dropboxusercontent.com/u/49022/diagnostictests/index.html
#' @export
#' @examples 
#' # Data.
#' dataset <- list(pop.size = 91, positives = 1)
#' 
#' # Initial conditions for chains.
#' inits <- list(list(true.prev = 0.05, se = 0.8, sp = 0.9),
#'               list(true.prev = 0.02, se = 0.3, sp = 0.7),
#'               list(true.prev = 0.09, se = 0.1, sp = 0.5))
#' 
#' # Priors.
#' priors <- c(true.prev.a = 1, true.prev.b = 1,
#'             se.a = 6.28, se.b = 13.32, sp.a = 212.12, sp.b = 3.13)
#' 
#' # Prevalence estimate.
#' prev.est <- OneTestOnePopBM(dataset = dataset, inits = inits, n.iter = 3e3,
#'                             priors = priors, pars = 'true.prev')
#' 
#' summary(prev.est)
#' 
#' # Diagnostic plots.
#' library(coda); library(ggmcmc)
#' gelman.diag(prev.est)
#' gelman.plot(prev.est)
#' gg.res <- ggs(prev.est)
#' ggs_traceplot(gg.res)
#' ggs_density(gg.res)
#' ggs_histogram(gg.res, bins = 100)
#' ggs_compare_partial(gg.res)
#' ggs_running(gg.res)
#' ggs_autocorrelation(gg.res)
OneTestOnePopBM <- function(dataset, inits, priors, pars, n.iter = 1e4,  n.chains = 3, burn.in = 1e3, thin = 1) {
    model <- 
        paste(c('model {',
                'positives ~ dbin(app.prev, pop.size)',
                'app.prev <- true.prev * se + (1 - true.prev) * (1 - sp)',
                'ppv <- true.prev * se / app.prev',
                'npv <- (1 - true.prev) * sp / (1 - app.prev)',
                paste0('se ~ dbeta(', priors[["se.a"]], ', ',
                       priors[["se.b"]], ')'),
                paste0('sp ~ dbeta(', priors[["sp.a"]], ', ',
                       priors[["sp.b"]], ')'),
                paste0('true.prev ~ dbeta(', priors[["true.prev.a"]], ', ',
                       priors[["true.prev.b"]], ')'),
                '}'),
              collapse = '\n')  
    writeLines(model, 'model.txt')
    jags.mod <- jags.model(file = 'model.txt',
                           data = dataset,
                           inits = inits, 
                           n.chains = n.chains)
    update(jags.mod, n.iter = burn.in)
    coda.samples(jags.mod,
                 variable.names = pars,
                 n.iter = n.iter,
                 thin = thin)
}