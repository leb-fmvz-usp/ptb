### http://blog.revolutionanalytics.com/2015/10/parameters-and-percentiles-the-gamma-distribution.html
# Parameters and percentiles (the gamma distribution)
# by Andrie de Vries
# 
# In one of John D. Cooke's blog posts of 2010 (Parameters and Percentiles), he poses the following problem:
# 
# The doctor says 10% of patients respond within 30 days of treatment and 80% respond within 90 days of treatment. Now go turn that into a probability distribution. That’s a common task in Bayesian statistics, capturing expert opinion in a mathematical form to create a prior distribution.
# 
# John then discusses how this level of information is highly valuable in statistical inference. The reason is that quite often this is the kind of information you might be able to elicit from a domain expert. It then is up to you as the statistician / data scientist to use this information.  In bayesian statistics, for example, you can use this information to construct a bayesian prior distribution. In particular, he demonstrates how this expectation can be modeled with a gamma distribution and shows how to solve the problem analytically.
# 
# In this post I demonstrate how to solve the problem using the non-linear least squares solver in R, using the nls() function.
# 
# But first, take a look at some of the properties of the gamma distribution. The gamma is a general family of distributions. Both the exponential and the chi-squared distributions are special cases of the gamma.
# 
# The gamma distribution takes two arguments. The first defines the shape. If shape is close to zero, the gamma is very similar to the exponential. If shape is large, then the gamma is similar to the chi-squared distribution.
# 
# To create the plots, you can use the function curve() to do the actual plotting, and dgamma() to compute the gamma density distribution. In this grid of plots, the shape parameter varies horisontally (from 1 on the left to 6 on the right). At the same time, the scale parameter varies vertically (from 0.1 at the top to 1.0 at the bottom).
# 
# 
# Next, you can use the function nls() to solve the problem as posed by John Cooke. The nls() function takes a loss function as an argument. This loss function is the function to be minimised by the solver. In the posed problem, you can compute the loss function as the difference between a hypothetical gamma distribution, calculated by qgamma() and the expected values posed by the problem.
# 
# The nls() solver is sensitive to the starting conditions, but easily finds a solution:







# Objective
#
# Fit a gamma distribution knowing that:
# - 20% fall below 15 days
# - 80% fall below 60 days
# Inspired by http://www.johndcook.com/blog/2010/01/31/parameters-from-percentiles/

x <- c(0.2, 0.8)
y <- c(15, 60)

# ?dgamma
# ?qgamma

# plot the gamma curve for a given shape and rate argument ----------------

# draws vertical lines at p
plotGamma <- function(shape=2, rate=0.5, to=0.99, p=c(0.1, 0.9), cex=1, ...){
  to <- qgamma(p=to, shape=shape, rate=rate)
  curve(dgamma(x, shape, rate), from=0, to=to, n=500, type="l", 
        main=sprintf("gamma(x, shape=%1.2f, rate=%1.2f)", shape, rate),
        bty="n", xaxs="i", yaxs="i", col="blue", xlab="", ylab="", 
        las=1, lwd=2, cex=cex, cex.axis=cex, cex.main=cex, ...)
  gx <- qgamma(p=p,  shape=shape, rate=rate)
  gy <- dgamma(x=gx, shape=shape, rate=rate)
  for(i in seq_along(p)) { lines(x=rep(gx[i], 2), y=c(0, gy[i]), col="blue") }
  for(i in seq_along(p)) { text(x=gx[i], 0, p[i], adj=c(1.1, -0.2), cex=cex) }
}

# plot several gamma curves -----------------------------------------------

# note the shape argument determines the overall shape
# the scale parameter only affects the scale values
oldpar <- par(mfrow = c(2, 3), mai=rep(0.5, 4))
plotGamma(1, 0.1)
plotGamma(2, 0.1)
plotGamma(6, 0.1)
plotGamma(1, 1)
plotGamma(2, 1)
plotGamma(6, 1)
par(oldpar)


# formulate a gamma error function for non-linear least squares -----------

errorGamma <- function(p=c(10, 3), quantiles, exp){
  gx <- qgamma(p=quantiles, shape=p[1], rate=p[2])
  sum((gx-exp)^2)
}

p <- c(0.1, 0.8)
solution <- nlm(f=errorGamma, p=c(2, 1), quantiles=p, exp=c(30, 90))

(shape <- solution$estimate[1])
(rate  <- solution$estimate[2])
plotGamma(shape, rate, p=p)
