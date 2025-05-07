# x<-X
prmtrs.init<-function (x, G, cluster, start = 1)
{
  n <- nrow(x)
  p <- ncol(x)
  a <- numeric(G)
  b <- matrix(NA, n, G)
  init.prmtrs <- list()
  wt <- matrix(0, n, G)
  init.prmtrs <- list(alpha = matrix(NA, p,G),
                      talpha = matrix(NA,G,p),
                      sig = array(NA, dim = c(p, p, G)),
                      inv.sig = array(NA,dim = c(p, p, G)),
                      log.det.sig = numeric(G),
                      mu = matrix(NA, p,G),
                      p = ncol(x),
                      n = nrow(x))
  for (g in 1:G) {
    wt[which(cluster == g), g] <- 1
    temp <- cov.wt(x = x, wt = wt[, g], center = TRUE, method = "ML")
    init.prmtrs$alpha[, g] <- rep(0, p)
    init.prmtrs$talpha[g, ] <- rep(0, p)
    init.prmtrs$mu[,g ] <- temp$center
    init.prmtrs$sig[, , g] <- temp$cov
    init.prmtrs$inv.sig[, , g] <- solve(init.prmtrs$sig[,, g])

    init.prmtrs$log.det.sig[g] <- log(det(matrix(init.prmtrs$sig[, , g],p,p)))
  }
  init.prmtrs$n.g <- apply(wt, 2, sum)
  init.prmtrs$pi.g <- apply(wt, 2, mean)
  init.prmtrs$n <- n
  init.prmtrs$p <- p
  init.prmtrs$nu <- (2 - p)/2
  return(init.prmtrs)
}


e1e2f <- function(p,x,mu,alpha,inv.sig){
  a <- 2 + c(alpha%*%inv.sig%*%alpha)
  mu <- as.vector(mu)
  b <- apply(x,1,function(v){(v-mu)%*%inv.sig%*%(v-mu)})
  # b<-apply(x, 1, function(v) {
  #   mu_i <- mu[1 , ]  # Select the correct column/component (adjust as needed)
  #   (v - mu_i) %*% inv.sig %*% (v - mu_i)
  # })
  t1 <- exp((log(a)+log(b))/2)
  t2 <- log(besselK(x=t1,nu=(2-p)/2 + 1,expon.scaled=TRUE)) - t1
  t3 <- log(besselK(x=t1,nu=(2-p)/2, expon.scaled=TRUE)) - t1
  t4 <- t2 - t3
  t5 <- c(log(b)-log(a))/2
  t6 <- sign((2-p)/2)*exp(log(2)+log(abs((2-p)/2)) - log(b))
  E1 <- exp(t5+t4)
  E2  <- exp(t4-t5) - t6
  out <- list();
  out$E1 <- E1
  out$E2 <- E2
  return(out)
}
# INF



e1e2fY <- function(p,x,mu,alpha,inv.sig){
  a <- 2 + c(alpha%*%inv.sig%*%alpha)
  mu <- as.matrix(mu)
  # b <- apply(x,1,function(v){(v-mu)%*%inv.sig%*%(v-mu)})
  b <- diag((x-as.vector(mu))%*%inv.sig%*%t(x-as.vector(mu)))
  # b <- ?mahalanobis(x = x, center = mu, cov = inv.sig, inverted = TRUE)

  # b<-apply(x, 1, function(v) {
  #   mu_i <- mu[1 , ]  # Select the correct column/component (adjust as needed)
  #   (v - mu_i) %*% inv.sig %*% (v - mu_i)
  # })
  t1 <- exp((log(a)+log(b))/2)
  t2 <- log(besselK(x=t1,nu=(2-p)/2 + 1,expon.scaled=TRUE)) - t1
  t3 <- log(besselK(x=t1,nu=(2-p)/2, expon.scaled=TRUE)) - t1
  t4 <- t2 - t3
  t5 <- c(log(b)-log(a))/2
  t6 <- sign((2-p)/2)*exp(log(2)+log(abs((2-p)/2)) - log(b))
  E1 <- exp(t5+t4)
  E2  <- exp(t4-t5) - t6
  out <- list();
  out$E1 <- E1
  out$E2 <- E2
  return(out)
}

check.val <- function(x,mu){
  v1 <- apply(x,1,function(v){ifelse(sqrt(sum((v-mu)^2))>1e-10 , 0, 1)})
  return(sum(v1))
}



check.valY <- function(x, mu) {
  # Ensure mu is a matrix
  v1 <- ifelse(sqrt(rowSums((x-mu)^2))>1e-6,0,1)
  return(v1)
}


dSALY <- function(X,mu,alpha,sigma){
  # print(sigma)
  X<-as.matrix(X)
  mu<- as.matrix(mu)
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
  ty.t <- X-as.vector(mu);
  inv.sig <- solve(sig, tol=1e-20)
  t1.num <- sweep(ty.t%*%inv.sig%*%talpha,2,log(2),"+");
  t1.den <- (p/2)*log(2*pi)+0.5*log.det.sig;
  t1 <- sweep(t1.num,2,t1.den,"-");
  maha <- diag((ty.t)%*%inv.sig%*%t(ty.t))
  # maha <- mahalanobis(x=X,center=mu,cov=inv.sig,inverted=TRUE);
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


clas <- function(groups,G){

  n <- length(groups)
  z <- array(0,c(n,G),dimnames=list(1:n,paste("comp",1:G,sep="")))
  for(i in 1:n)
    z[i,groups[i]] <- 1
  return(z)

}

aitkens<-function (lval, i, eps)
{
  if (i > 2) {
    aa <- (lval[i] - lval[i - 1])/(lval[i - 1] - lval[i -
                                                        2])
    inf.l <- lval[i - 1] + (1/(1 - aa)) * (lval[i] - lval[i -
                                                            1])
    return(inf.l)
  }
  else {
    inf.l <- lval[i] + eps + 1
    return(inf.l)
  }
}
