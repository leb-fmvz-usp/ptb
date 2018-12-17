library(R2OpenBUGS)

## One test, one population binomial model -------------------------------------

## Terms
# ------------------------------------------------
# Term                      | ptb       | OpenBUGS
# ------------------------------------------------
# True prevalence           | true_prev | pi
# Apparent prevalence       | app_prev  | p
# Sensitivity               | se        | se
# Specificity               | sp        | sp
# Positive predicitve value | ppv       | ppv
# Negative predicitve value | npv       | npv
# ------------------------------------------------

## ptb
dataset <- list(pop_size = 91, positives = 1)
priors <- c(true_prev_a = 1.8, true_prev_b = 26.74,
            se_a = 22.5, se_b = 10.22,
            sp_a = 88.28, sp_b = 1.882)
pars <- c("true_prev", "app_prev", "ppv", "npv")
prev_est <- OneTestOnePopBM(dataset = dataset, n_iter = 3e3,
                            priors = priors,
                            pars = pars)

## OpenBUGS
dataset_ob <- c(n = 91, y = 1)
pars_ob <- c("pi", "p", "ppv", "npv")
prev_est_ob <- bugs(dataset_ob,
                    inits = NULL,
                    model.file = "one_test_one_pop_bm.txt",
                    parameters = pars_ob,
                    n.chains = 3,
                    n.burnin = 1000,
                    n.iter = 3e3)

summary(prev_est)
print(prev_est_ob, digits.summary = 5)

## One test, one population binomial mixture model -----------------------------

## Terms
# ---------------------------------------------------------------------
# Term                                        | ptb            | OpenBUGS
# ---------------------------------------------------------------------
# True prevalence                             | true_prev      | pi
# True prevalence within positive herds       | true_prev_wph  | pistar
# Herd prevalence                             | prev_h         | z
# ---------------------------------------------------------------------

## ptb
dataset <- list(pop_size = 123, positives = 3)
priors <- list(true_prev_wph_a = 1, true_prev_wph_b = 1,
               se_a = 103, se_b = 118,
               sp_a = 275, sp_b = 40,
               prev_h = 0.9)
prev_est <- OneTestOnePopBMM(dataset = dataset, priors = priors, n_iter = 3e3,
                             pars = c('true_prev', 'true_prev_wph', 'prev_h'))

## OpenBUGS
dataset_ob <- c(n = 123, y = 3, inf = .9)
pars_ob <- c("pi", "pistar", "z")
prev_est_ob <- bugs(dataset_ob,
                    inits = NULL,
                    model.file = "one_test_one_pop_bmm.txt",
                    parameters = pars_ob,
                    n.chains = 3,
                    n.burnin = 1000,
                    n.iter = 3e3)

summary(prev_est)
print(prev_est_ob, digits.summary = 5)

## One test, multiple populations binomial mixture model -----------------------

## Terms
# ---------------------------------------------------------------------
# Term                                        | ptb           | OpenBUGS
# ---------------------------------------------------------------------
# Population size per herd                    | pop_size      | n
# Number of positive animals per herd         | positives     | y
# Probability of the herd being positive      | prev_h        | tau
# True prevalence within positive herds       | true_prev_wph | pi
# ---------------------------------------------------------------------

## ptb
q1 <- ElicitBeta(.3, .99, .5, summary = T, quantiles = c(.95))
curve(dbeta(x, q1[[1]], q1[[2]]), 0, 1)
dataset <- list(n = 29,
                pop_size = rep(60, 29),
                positives = c(2, 1, 2, 2, 3, 6, 0, 6, 3, 13, 2, 3, 1, 7,
                              2, 2, 0, 4, 1, 2, 6, 1, 4, 0, 5, 4, 2, 0, 13))
priors <- list(true_prev_wph_mean_a = 3.283, true_prev_wph_mean_b = 17.744,
               se_a = 58.8, se_b = 174.5, sp_a = 272.4, sp_b = 6.5,
               prev_h_a = 4.8, prev_h_b = 3.6,
               true_prev_wph_var_a = 2, true_prev_wph_var_b = 0.1)
pars <- c('pred_true_prev_h',
          'pred_true_prev_wph',
          'pred_inf_h',
          'prev_h')
prev_est <- OneTestMultPopBMM(dataset = dataset, priors = priors, n_iter = 5e3,
                              pars = pars)

## OpenBUGS
dataset_ob <- list(n=60,
                   y = c(2, 1, 2, 2, 3, 6, 0, 6, 3, 13, 2, 3, 1, 7,
                         2, 2, 0, 4, 1, 2, 6, 1, 4, 0, 5, 4, 2, 0, 13))
pars_ob <- c("pi30", "pistar30", "Z30", "tau", "a1", "a2", "a3", "b1")
#pars_ob <- c("Se", "Sp", "tau")
prev_est_ob <- bugs(dataset_ob,
                    inits = NULL,
                    model.file = "one_test_mult_pop_bmm2.txt",
                    parameters = pars_ob,
                    n.chains = 3,
                    n.burnin = 1000,
                    n.iter = 5e3)

summary(prev_est)
print(prev_est_ob, digits.summary = 5)
str(prev_est)
mean(sapply(prev_est, function(x) sum(x[, 2] < .05) / 5000))
mean(sapply(prev_est, function(x) sum(x[, 2] < .5) / 5000))
mean(sapply(prev_est, function(x) sum(x[, 2] == 0) / 5000))
mean(sapply(prev_est, function(x) sum(x[, 4] > .5) / 5000))

q1 <- function(model = NULL, lessthan = NULL, greaterthan = NULL, equal = NULL, par = NULL) {
    if (!is.null(lessthan)) {
        return(mean(sapply(model, function(x) sum(x[, par] < lessthan) / 5000)))
    }
    if (!is.null(equal)) {
        return(mean(sapply(model, function(x) sum(x[, par] == equal) / 5000)))
    }
    if (!is.null(greaterthan)) {
        return(mean(sapply(model, function(x) sum(x[, par] > greaterthan) / 5000)))
    }
}
q1(prev_est, lessthan = .05, par = "pred_true_prev_h")
q1(prev_est, lessthan = .5, par = "pred_true_prev_h")
q1(prev_est, equal = 0, par = "pred_true_prev_h")
q1(prev_est, greaterthan = .5, par = "prev_h")
