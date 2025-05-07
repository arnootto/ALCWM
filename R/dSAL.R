#' shifted asymmetric Laplace density.
#'
#' Density function of the shifted asymmetri Laplace distribution.
#'
#' @param x A n by p matrix where each row corresponds a p-dimensional observation.
#' @param mu A vector specifying the mean vector.
#' @param alpha A vector specifying the direction of each skewness in each variable.
#' @param sigma A matrix specifying the covariance matrix of the variables.
#'
#' @return A vector of length n that gives the value of the probability density function for each obervation in the matrix x and the specified parameter values.
#'
#' @import stats
#'
#' @export
dSAL <- function(X,mu,alpha,sigma){
  # print(sigma)
  p <- ncol(X)
  sig=sigma
  x <- as.matrix(X)
  talpha <- as.matrix(alpha,nrow=1,ncol=p);
  alpha <- t(talpha);
  n <- nrow(x);
  nu <- (2-p)/2;
  if(p>1){ log.det.sig <- log(det(sig))
  }else{log.det.sig=log(sig)}
  # if(log.det.sig == "NaN") stop('dsal Degenerate Solution Reached - The Determinate of this matrix is Zero');
  ty.t <- sweep(x,2,mu,"-");
  inv.sig <- solve(sig, tol=1e-20)
  t1.num <- sweep(ty.t%*%inv.sig%*%talpha,2,log(2),"+");
  t1.den <- (p/2)*log(2*pi)+0.5*log.det.sig;
  t1 <- sweep(t1.num,2,t1.den,"-");
  # diag((x-mu)%*%inv.sig%*%t(x-mu))
  maha <- mahalanobis(x=x,center=mu,cov=inv.sig,inverted=TRUE);
  ata <- 2 + alpha%*%inv.sig%*%talpha;
  l.t2.num <- log(maha)
  l.t2.den <- log(ata)
  l.t2.den <- rep(l.t2.den,n)
  t2 <- 0.5*nu*(l.t2.num-l.t2.den)
  u <- exp( 0.5*(l.t2.den + l.t2.num) )
  t3 <- log(besselK(u,nu,expon.scaled=TRUE)) - u
  val1 <- t1+t2+t3
  val <- exp(val1)
  return(c(val))
}
