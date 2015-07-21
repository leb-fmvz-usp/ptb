#' @title 1 test and 1 population binomial mixture model
#' @description 1 test and 1 population binomial mixture model to estimate prevalence and diagnostic-test-related meassures.
#' @param dataset \code{\link{list}} with the population size and the number of positives. The names of these values must be "pop.size" and "positives" respectively.
#' @param inits \code{\link{list}} with initial conditions for chains. \code{inits} must define initial values for the true prevalence within infected herds, the sensitivity, the specificity and the herd prevalence. The names of these values must be "true.prev.wph", "se", "sp" and "prev.h" respectively.
#' @param priors \code{\link{vector}} with the parameters a and b (Beta distribution) for the true prevalence within infected herds, the sensitivity and the specificity; and with the parameter p (Bernoulli distribution) for the herd prevalence. The names of these values must be: "true.prev.wph.a", "true.prev.wph.b", "se.a", "se.b", "sp.a", "sp.b" and "prev.h".
#' @param pars character vector giving the names of parameters to be monitored. It is passed to the \code{variable.names} argument of the \code{\link{coda.samples}} function. In addition to the parameters specified in \code{inits}, the apparent prevalence, the true prevalence, and the positive and negative predictive values can be monitored if specified in \code{pars} as "app.prev", "true.prev", "ppv" and "npv", respectively.
#' @param n.chains the number of parallel chains for the model. It is passed to the \code{n.chains} argument of the \code{\link{jags.model}} function.
#' @param burn.in the number of iteration to be discarded. It is passed to the \code{n.iter} argument of the \code{\link{update.jags}} function.
#' @param thin thinning interval for monitors. It is passed to the \code{thin} argument of the \code{\link{coda.samples}} function.
#' @param n.iter number of iterations to monitor. It is passed to the \code{n.iter} argument of the \code{\link{coda.samples}} function.
#' @return A \code{\link{list}} of class \code{mcmc.list}.
#' @references https://dl.dropboxusercontent.com/u/49022/diagnostictests/index.html
#' @export
#' @examples 
#' # Data (initial values for chains automatically generated).
#' dataset <- list(pop.size = 91, positives = 1)
#' 
#' # Priors.
#' priors <- list(true.prev.wph.a = 1.8, true.prev.wph.b = 26.74, se.a = 6.28,
#'                se.b = 13.32, sp.a = 212.12, sp.b = 3.13, prev.h = 0.1)
#'                
#' 
#' # Prevalence estimates.
#' prev.est <- OneTestOnePopBMM(dataset = dataset, priors = priors, n.iter = 3e3,
#'                              pars = c('true.prev', 'true.prev.wph', 'prev.h'))
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
OneTestOnePopBMM <- function(dataset, inits, priors, pars, n.iter = 1e4, n.chains = 3, burn.in = 1e3, thin = 1) {
    model <- 
        paste(c('model {',
                'positives ~ dbin(app.prev, pop.size)',
                'app.prev <-true.prev * se + (1 - true.prev) * (1 - sp)',
                'ppv <- true.prev * se / app.prev',
                'npv <- (1 - true.prev) * sp / (1 - app.prev)',
                paste0('prev.h ~ dbern(', priors[["prev.h"]], ')'),
                paste0('se ~ dbeta(', priors[["se.a"]], ', ',
                       priors[["se.b"]], ')'),
                paste0('sp ~ dbeta(', priors[["sp.a"]], ', ',
                       priors[["sp.b"]], ')'),
                paste0('true.prev.wph ~ dbeta(',
                       priors[['true.prev.wph.a']],
                       ', ', priors[['true.prev.wph.b']], ')'),
                'true.prev <- prev.h * true.prev.wph',
                '}'),
              collapse = '\n')
    writeLines(model, 'model.txt')
    jags.mod <- jags.model(file = 'model.txt',
                          data = dataset,
                          inits = inits, 
                          n.chains = 3)
    update(jags.mod, n.iter = burn.in)
    coda.samples(jags.mod,
                 variable.names = pars,
                 n.iter = n.iter,
                 thin = thin)
}