#' @title Two dependent tests, two populations and no gold standard
#' @description 2 tests and 1 population binomial model to estimate prevalence and diagnostic test related meassures.
#' @param dataset \code{\link{list}} with the population size and test results. The names of these values must be "pop_size" and "t_res" respectively. \code{t_res} is a named vector with the number of animals thesting postive in both tests "t1p_t2p", testing positive in the first and negative in the second "t1p_t2n", testing negative in the second and positive in the second "t1n_t2p", and testing positive in both "t1n_t2n".
#' @param inits \code{\link{list}} with initial conditions for chains. \code{inits} must define initial values for pi, se_test1, sp_test1, se_test2 and sp_test2.
#' @param priors \code{\link{vector}} with the parameters a and b (Beta distribution) for pi, se_test1, sp_test1, se_test2 and sp_test2..
#' @param pars character vector giving the names of parameters to be monitored. It is passed to the \code{variable.names} argument of the \code{\link{coda.samples}} function. In addition to the parameters specified in \code{inits}, ... monitored if specified in \code{pars} as "app_prev", "ppv" and "npv", respectively.
#' @param n_chains the number of parallel chains for the model. It is passed to the \code{n.chains} argument of the \code{\link{jags.model}} function.
#' @param burn_in the number of iteration to be discarded. It is passed to the \code{n.iter} argument of the \code{\link{update.jags}} function.
#' @param thin thinning interval for monitors. It is passed to the \code{thin} argument of the \code{\link{coda.samples}} function.
#' @param n_iter number of iterations to monitor. It is passed to the \code{n.iter} argument of the \code{\link{coda.samples}} function.
#' @details This function creates a text file with the model and it is saved in the working directory.
#' @return A \code{\link{list}} of class \code{mcmc.list}.
#' @references https://cadms.vetmed.ucdavis.edu/diagnostic/software
#' @export
#' @examples 
#' # Dataset
#' dataset <- list(pop_size = 214,
#'                 t_res = c(t1p_t2p = 121, t1p_t2_n = 6,
#'                           t1n_t2_p = 16, t1n_t2n = 71))
#' 
#' # Priors
#' priors <- c(pi_a = 13.322, pi_b = 6.281,
#'             se_test1_a = 9.628, se_test1_b = 3.876,
#'             sp_test1_a = 15.034, sp_test1_b = 2.559,
#'             se_test2_a = 9.628, se_test2_b = 3.876,
#'             sp_test2_a = 15.034, sp_test2_b = 2.559)
#' 
#' # Estimates
#' est <- TwoDepTestsOnePopNGS(dataset = dataset, n_iter = 3e3,
#'                                  priors = priors, pars = c("se_test1", "se_test2"))
#' 
#' summary(est)
#'
#' # Diagnostic plots.
#' library(coda); library(ggmcmc)
#' gelman.diag(est)
#' gelman.plot(est)
#' gg_res <- ggs(est)
#' ggs_traceplot(gg_res)
#' ggs_density(gg_res)
#' ggs_histogram(gg_res, bins = 100)
#' ggs_compare_partial(gg_res)
#' ggs_running(gg_res)
#' ggs_autocorrelation(gg_res)
#' 
TwoDepTestsOnePopNGS <- function(dataset, inits, priors, pars, n_iter = 1e4,  n_chains = 3, burn_in = 1e3, thin = 1) {
  model <- 
    paste(c("model{",
            "t_res[1:4] ~ dmulti(p[1:4], pop_size)",
            "p[1] <- pi*(se_test1*se_test2+covDp) + (1-pi)*((1-sp_test1)*(1-sp_test2)+covDn)",
            "p[2] <- pi*(se_test1*(1-se_test2)-covDp) + (1-pi)*((1-sp_test1)*sp_test2-covDn)",
            "p[3] <- pi*((1-se_test1)*se_test2-covDp) + (1-pi)*(sp_test1*(1-sp_test2)-covDn)",
            "p[4] <- pi*((1-se_test1)*(1-se_test2)+covDp) + (1-pi)*(sp_test1*sp_test2+covDn)",
            "ls <- (se_test1-1)*(1-se_test2)",
            "us <- min(se_test1,se_test2) - se_test1*se_test2",
            "lc <- (sp_test1-1)*(1-sp_test2)",
            "uc <- min(sp_test1,sp_test2) - sp_test1*sp_test2",
            paste0('pi ~ dbeta(', priors[["pi_a"]], ', ',
                   priors[["pi_b"]], ')'),
            paste0('se_test1 ~ dbeta(', priors[["se_test1_a"]], ', ',
                   priors[["se_test1_b"]], ')'),
            paste0('sp_test1 ~ dbeta(', priors[["sp_test1_a"]], ', ',
                   priors[["sp_test1_b"]], ')'),
            paste0('se_test2 ~ dbeta(', priors[["se_test2_a"]], ', ',
                   priors[["se_test2_b"]], ')'),
            paste0('sp_test2 ~ dbeta(', priors[["sp_test2_a"]], ', ',
                   priors[["sp_test2_b"]], ')'),
            "covDn ~ dunif(lc, uc)",
            "covDp ~ dunif(ls, us)",
            "rhoD <- covDp / sqrt(se_test1*(1-se_test1)*se_test2*(1-se_test2))",
            "rhoDc <- covDn / sqrt(sp_test1*(1-sp_test1)*sp_test2*(1-sp_test2))",
            "}"),
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
