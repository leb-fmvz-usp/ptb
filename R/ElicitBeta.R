#' @title Parameterization of a Beta distribution
#' @description Parameterization of a Beta distribution from elicited mode and quantiles.
#' @param mode scalar representing the elicited mode of the variable of interest. Must be a number between 0 and 1.
#' @param maximum scalar representing the elicited maximum value of the variable of interest. Must be a number between 0 and 1.
#' @param minimum scalar representing the elicited minimum value of the variable of interest. Must be a number between 0 and 1.
#' @param confidence scalar representing the the level of confidence on the definition of \code{maximum} or \code{minimum}. Must be a number between 0 and 1.
#' @param summary logical. If FALSE (default), only the parameters a (shape1) and b (shape2) are returned. If TRUE, the mean the variance and the quantiles are also returned.
#' @param quantiles \code{\link{numeric}} \code{\link{vector}} with the quantiles to be computed if \code{summary == TRUE}. The quantiles must be numbers between 0 and 1.
#' @details The parameters a and b are calculated using the algorithm implemented in the \code{epi.betabuster} function of the epiR package.
#' @return If \code{summary == FALSE}, a \code{\link{numeric}} \code{\link{vector}} with the prameters a and b. If \code{summary == TRUE}, a \code{\link{list}} with the parameters a and b, as well as summary statistics of the distribution.
#' @export
#' @examples
#' ElicitBeta(mode = 0.3, maximum = 0.5, confidence = 0.95, summary = FALSE)
#' ElicitBeta(mode = 0.99, minimum = 0.97, confidence = 0.95, summary = TRUE)
ElicitBeta <- function (mode, confidence, maximum = NULL, minimum = NULL, summary = FALSE, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
    if (mode < 0 | mode > 1) {
        stop('mode must be a number between 0 and 1.')
    }
    if (mode < 0 | mode > 1) {
        stop('confidence must be a number between 0 and 1.')
    }
    if (!is.null(minimum) & !is.null(maximum)) {
        stop('Define maximum or minimum, not both.')
    }
    if (is.null(minimum) & is.null(maximum)) {
        stop('Define maximum or minimum.')
    }
    a <- seq(1, 1e4, by = 0.01)
    b <- 2 - a + (a - 1) / mode
    if (!is.null(maximum)) {
        if (maximum < 0 | maximum > 1) {
            stop('maximum must be a number between 0 and 1.')
        }
        x <- maximum
        dist <- pbeta(q = x, shape1 = a, shape2 = b)
        idx <- which((abs(dist - confidence)) == min(abs(dist - confidence)))
    }
    if (!is.null(minimum)) {
        if (minimum < 0 | minimum > 1) {
            stop('maximum must be a number between 0 and 1.')
        }
        x <- minimum
        dist <- pbeta(q = x, shape1 = a, shape2 = b)
        idx <- which((abs(dist - (1 - confidence))) ==
                         min(abs(dist - (1 - confidence))))
    }
    if (!summary) {
        res <- c('a (shape1)' = a[idx], 'b (shape2)' = b[idx])
        class(res) <- 'ElicitBeta'
        return(res)
    }
    if (summary) {
        if (min(quantiles) < 0 | max(quantiles) > 1) {
            stop('quantiles must be numbers between 0 and 1.')
        }
        a = a[idx]
        b = b[idx]
        beta.mean <- a / (a + b)
        beta.var <- a * b / (((a + b) ^ 2) * (a + b + 1))
        beta.qts <- qbeta(quantiles,shape1 = a, shape2 = b)
        names(beta.qts) <- quantiles
        res <- list('a (shape1)' = a, 'b (shape2)' = b, Mean = beta.mean,
                    Variance = beta.var, Quantiles = beta.qts)
        class(res) <- 'ElicitBeta'
        return(res)
    }
}
