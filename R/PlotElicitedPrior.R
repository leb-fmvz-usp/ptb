#' Plot elicited priors.
#' @description Plots a Beta distribution parameterized with the \code{\link{ElicitBeta}} function.
#' @param elicited.prior output form the \code{\link{ElicitBeta}} function.
#' @param ylab a title for the y axis.
#' @param ... further arguments passed to the \code{\link{curve}} fuction.
#' @return Density plot.
#' @export
#' @examples 
#' elicited.prior <- ElicitBeta(mode = 0.3, maximum = 0.5, confidence = 0.95)
#' PlotElicitedPrior(elicited.prior, col = 'red')
PlotElicitedPrior <- function(elicited.prior, ylab = 'Density', ...) {
    if (class(elicited.prior) != 'ElicitBeta') {
        stop('The class of elicited.prior must be ElicitBeta')
    }
    x <- c()
    curve((dbeta(x, as.numeric(elicited.prior[1]),
                 as.numeric(elicited.prior[2]))), ylab = ylab, ...)
}
