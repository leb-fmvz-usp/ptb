#' @title 1 test and multiple populations binomial mixture model
#' @description 1 test and multiple population binomial mixture model to estimate herd prevalence and prevalence within positive herds.
#' @param dataset \code{\link{list}} with the population size and the number of positives. The names of these values must be "pop_size" and "positives" respectively.
#' @param inits \code{\link{list}} with initial conditions for chains. \code{inits} must define initial values for the true prevalence within infected herds, the sensitivity, the specificity and the herd prevalence. The names of these values must be "true_prev_wph", "se", "sp" and "prev_h" respectively.
#' @param priors \code{\link{vector}} with the parameters a and b (Beta distribution) for the true prevalence within infected herds, the sensitivity and the specificity; and with the parameter p (Bernoulli distribution) for the herd prevalence. The names of these values must be: "true_prev_wph_mean_a", "true_prev_wph_mean_b", "se_a", "se_b", "sp_a", "sp_b" "prev_h_a", "prev_h_b", true_prev_wph_var_a", "true_prev_wph_var_b".
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
OneTestMultPopBMM <- function(dataset, inits, priors, pars, n_iter = 1e4, n_chains = 3, burn_in = 1e3, thin = 1) {
    model <-
        paste(c('model {',
                'for (i in 1:n) {',
                    paste(c('app_prev[i] <-true_prev_h[i] * se + (1 - true_prev_h[i]) * (1 - sp)',
                            'positives[i] ~ dbin(app_prev[i], pop_size[i])',
                            paste0('inf_h[i] ~ dbern(prev_h)'),
                            'true_prev_wph[i] ~ dbeta(true_prev_wph_a, true_prev_wph_b)',
                            'true_prev_h[i] <- true_prev_wph[i] * inf_h[i]'),
                          collapse = "\n"),
                '}',
                paste0('se ~ dbeta(', priors[["se_a"]], ', ',
                       priors[["se_b"]], ')'),
                paste0('sp ~ dbeta(', priors[["sp_a"]], ', ',
                       priors[["sp_b"]], ')'),
                paste0('prev_h ~ dbeta(', priors[["prev_h_a"]], ', ',
                       priors[["prev_h_b"]], ')'),
                'true_prev_wph_a <- true_prev_wph_mean * true_prev_wph_var',
                'true_prev_wph_b <- true_prev_wph_var * (1 - true_prev_wph_mean)',
                paste0('true_prev_wph_mean ~ dbeta(',
                       priors[["true_prev_wph_mean_a"]], ', ',
                       priors[["true_prev_wph_mean_b"]], ')'),
                paste0('true_prev_wph_var ~ dgamma(',
                       priors[["true_prev_wph_var_a"]], ', ',
                       priors[["true_prev_wph_var_b"]], ')'),
                paste(c('pred_inf_h ~ dbern(prev_h)',
                        'pred_true_prev_wph ~ dbeta(true_prev_wph_a, true_prev_wph_b)',
                        'pred_true_prev_h <- pred_true_prev_wph * pred_inf_h'),
                      collapse = "\n"),
                # 'for (i in 1:sims) {',
                #     paste(c('pred_inf_h[i] ~ dbern(prev_h)',
                #             'pred_true_prev_wph[i] ~ dbeta(true_prev_wph_a, true_prev_wph_b)',
                #             'pred_true_prev_h[i] <- pred_true_prev_wph[i] * pred_inf_h[i]'),
                #           collapse = "\n"),
                # '}',
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