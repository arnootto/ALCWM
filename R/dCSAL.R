#' contaminated shifted asymmetric Laplace density.
#'
#' Density function of the contaminated shifted asymmetri Laplace distribution.
#'
#' @param x A n by p matrix where each row corresponds a p-dimensional observation.
#' @param mu A vector specifying the mean vector.
#' @param alpha A vector specifying the direction of each skewness in each variable.
#' @param sigma A matrix specifying the covariance matrix of the variables.
#' @param delta proportion of atypical values. 0 < delta < 1.
#' @param eta inflation parameter. eta > 1.
#'
#' @return A vector of length n that gives the value of the contaminated shifted asymmetric Laplace probability density function for each obervation in the matrix x and the specified parameter values.
#'
#' @import stats
#'
#' @export
dCSAL <- function(X,mu=rep(0,p),alpha=rep(0,p),sigma=diag(p),delta=0.9,eta=1.01) {
  dg=dSAL(X,alpha=alpha,sigma=sigma,mu=mu)  #dd=sweep(data%*%gam,2,mu%*%gam+sqrt(rho)*skew%*%gam,FUN="-")
  db=dSAL(X,alpha=(sqrt(eta)*alpha),sigma=(eta*sigma),mu=mu)
  den = delta*dg+(1-delta)*db
  return(den)
}
