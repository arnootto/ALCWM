#' Contaminated Shifted Asymmetric Laplace Cluster-Weighted Model
#'
#' Maximum likelihood fitting of the contaminated shifted asymmetric Laplace (SAL) cluster-weighted model (CWM) using an ECM algorithm.
#'
#' @param Y A matrix or vector containing the response variables.
#' @param X A matrix or vector containing the predictor variables.
#' @param G Number of clusters (components) to fit. Default is 2.
#' @param tol Convergence tolerance for the outer EM algorithm. Default is 0.01.
#' @param max.it Maximum number of iterations for the outer EM algorithm. Default is 2000.
#' @param salcwm.tol Convergence tolerance for the inner SALCWM algorithm. Default is 0.01.
#' @param salcwm.max.it Maximum number of iterations for the inner SALCWM algorithm. Default is 100.
#' @param etamax Upper bound for the inflation parameter (eta). Default is 100.
#' @param initialization Method for initializing the posterior probabilities. Options are
#'   \code{"mclust"} (default), \code{"kmeans"}, \code{"random.soft"}, \code{"random.hard"},
#'   or \code{"manual"}.
#' @param print.iter Logical. If \code{TRUE}, prints the iteration number during fitting.
#'   Default is \code{TRUE}.
#' @param start.z Optional matrix of initial posterior probabilities for
#'   \code{initialization = "manual"}. Default is \code{NULL}.
#'   @param mu.tol Convergence tolerance for the mean.
#'
#' @return A list containing the following components:
#'   \itemize{
#'     \item \code{Y}: The input response matrix.
#'     \item \code{X}: The input predictor matrix.
#'     \item \code{G}: Number of clusters.
#'     \item \code{n}: Number of observations.
#'     \item \code{salcwm}: Results from the inner SALCWM model fit.
#'     \item \code{npar}: Number of model parameters.
#'     \item \code{prior}: Estimated prior probabilities for each cluster.
#'     \item \code{muX}: Estimated mean vectors for \code{X} in each cluster.
#'     \item \code{SigmaX}: Estimated covariance matrices for \code{X} in each cluster.
#'     \item \code{alphaX}: Estimated skewness parameters for \code{X} in each cluster.
#'     \item \code{deltaX}: Estimated proportion of atypical values for \code{X} in each cluster.
#'     \item \code{etaX}: Estimated inflation parameters for \code{X} in each cluster.
#'     \item \code{beta}: Estimated regression coefficients for \code{Y} given \code{X} in each cluster.
#'     \item \code{muY}: Estimated mean vectors for \code{Y} in each cluster.
#'     \item \code{SigmaY}: Estimated covariance matrices for \code{Y} in each cluster.
#'     \item \code{alphaY}: Estimated skewness parameters for \code{Y} in each cluster.
#'     \item \code{deltaY}: Estimated proportion of atypical values for \code{Y} in each cluster.
#'     \item \code{etaY}: Estimated inflation parameters for \code{Y} in each cluster.
#'     \item \code{z}: Posterior probabilities for each observation in each cluster.
#'     \item \code{group}: Cluster assignments based on maximum posterior probability.
#'     \item \code{detection}: Indicators for detected atypical observations.
#'     \item \code{density}: Mixture density values for each observation.
#'     \item \code{iter.stop}: Number of iterations until convergence.
#'     \item \code{loglik}: Final log-likelihood value.
#'     \item \code{AIC}: Akaike Information Criterion.
#'     \item \code{BIC}: Bayesian Information Criterion.
#'   }
#'
#'
#' @import mclust
#'
#'#'@examples
#'# Fit the CSALCWM model with G = 2 components to the AIS dataset
#'library(sn)
#'data("ais")
#'est=CSALCWM(Y=cbind(ais$RCC,ais$WCC),X=cbind(ais$BMI,ais$SSF,ais$Bfat,ais$LBM),G=2, tol=1e-5,max.it=2000,initialization = "mclust")
#
#'
#' @export
CSALCWM <-function(Y,X,G=2,
                   tol=1e-2,
                   max.it=2000,
                   salcwm.tol=1e-2,
                   salcwm.max.it=100,
                   etamax=100,
                   initialization = "mclust",
                   print.iter=T,
                   start.z=NULL,
                   mu.tol=1e-5){
  if(is.vector(Y))
    Y <- matrix(Y,ncol=1)
  if(is.vector(X))
    X <- matrix(X,ncol=1)
  if(any(is.na(X)) | any(is.na(Y)))
    stop('No NAs allowed.')
  dY <- ncol(Y)
  dX <- ncol(X)
  x<-X
  n  <- nrow(X)

  z  <- array(0,c(n,G),dimnames=list(1:n,paste("comp.",1:G,sep="")))

  npar <- G-1 + G * (dX + dX * (dX + 1)/2 + dX+2) + G * (dY * (dX + 1) + dY * (dY + 1)/2 + dY+2)

  prior  <- numeric(G)

  # X Parameters definition

  salcwm <- SALCWM(Y=Y,X=X,G=G,tol=salcwm.tol,max.it = salcwm.max.it, initialization = initialization, print.iter=print.iter, start.z=start.z, mu.tol=mu.tol)
  prmtrs <- salcwm$prmtrs
  zig.hat <- salcwm$z
  salcwm$loglik

  Xstar <- cbind(rep(1,n),X)

  # Distribution

  dens  <- array(0,c(n,G),dimnames=list(1:n,paste("comp.",1:G,sep="")))
  densX  <- array(0,c(n,G),dimnames=list(1:n,paste("comp.",1:G,sep="")))
  densY  <- array(0,c(n,G),dimnames=list(1:n,paste("comp.",1:G,sep="")))

  E1.X <- E2.X <- E1.Y <- E2.Y <- matrix(NA,n,G)
  E1.tilde.X <- E2.tilde.X <- E1.tilde.Y <- E2.tilde.Y <- matrix(NA,n,G)
  Good.distX <- Bad.distX <- Good.distY <- Bad.distY <- matrix(NA,n,G)
  Group.distX <- Group.distY <- Mix.dist <- matrix(NA,n,G)

  deltaX <- rep(0.99,G)
  etaX   <- rep(1.01,G)

  deltaY <- rep(0.99,G)
  etaY   <- rep(1.01,G)

  vig <-matrix(deltaX,n,G)
  uig <-matrix(deltaY,n,G)

  loglik <- NULL
  loglikX<- NULL
  loglikY<-NULL
  i=1
  curr.muX <- matrix(NA,dX,G)
  curr.muY <- prmtrs$muY
  Xstar <- cbind(1,X)
  a1 <- 0
  a2 <- tol + runif(1, 0, 1)

  while (abs(a1 - a2) > tol) {
    nu <- prmtrs$nu
    pi.g <- prmtrs$pi.g

    for (g in 1:G){
      e1e2 <- e1e2f(p=dX,x=X,mu=prmtrs$X$mu[,g],alpha=prmtrs$X$alpha[,g],inv.sig=prmtrs$X$inv.sig[,,g])
      E1.X[,g] <- e1e2$E1
      E2.X[,g] <- e1e2$E2
      e1e2.tilde <- e1e2f(p=dX,x=X,mu=prmtrs$X$mu[,g],alpha=sqrt(etaX[g])*prmtrs$X$alpha[,g],inv.sig=(1/etaX[g])*prmtrs$X$inv.sig[,,g])
      E1.tilde.X[,g] <- e1e2.tilde$E1
      E2.tilde.X[,g] <- e1e2.tilde$E2

      e1e2 <- e1e2fY(p=dY,x=Y,mu=prmtrs$muY[,,g],alpha=prmtrs$Y$alpha[,g],inv.sig=prmtrs$Y$inv.sig[,,g])
      E1.Y[,g] <- e1e2$E1
      E2.Y[,g] <- e1e2$E2
      e1e2.tilde <- e1e2fY(p=dY,x=Y,mu=prmtrs$muY[,,g],alpha=sqrt(etaY[g])*prmtrs$Y$alpha[,g],inv.sig=(1/etaY[g])*prmtrs$Y$inv.sig[,,g])
      E1.tilde.Y[,g] <- e1e2.tilde$E1
      E2.tilde.Y[,g] <- e1e2.tilde$E2

    }

    for (g in 1:G){
      Good.distX[,g] <- deltaX[g] * dSAL(X=X,mu=prmtrs$X$mu[,g],alpha=prmtrs$X$alpha[,g],sigma = prmtrs$X$sig[,,g])
      Good.distX[,g] <- sapply(Good.distX[,g],function(v){ifelse(v<2.205e-308,2.205e-308,v)})
      Bad.distX[,g] <- dSAL(X=X,mu=prmtrs$X$mu[,g],alpha=sqrt(etaX[g])*prmtrs$X$alpha[,g],sigma=etaX[g]*prmtrs$X$sig[,,g])
      Bad.distX[,g] <- sapply(Bad.distX[,g],function(v){ifelse(v<2.205e-308,2.205e-308,v)})
      Group.distX[,g] <- Good.distX[,g]+(1-deltaX[g])*Bad.distX[,g]


      Good.distY[,g] <- deltaY[g] * dSALY(X=Y,mu=prmtrs$muY[,,g],alpha=prmtrs$Y$alpha[,g],sigma = prmtrs$Y$sig[,,g])
      Good.distY[,g] <- sapply(Good.distY[,g],function(v){ifelse(v<2.205e-308,2.205e-308,v)})
      Bad.distY[,g] <- dSALY(X=Y,mu=prmtrs$muY[,,g],alpha=sqrt(etaY[g])*prmtrs$Y$alpha[,g],sigma=etaY[g]*prmtrs$Y$sig[,,g])
      Bad.distY[,g] <- sapply(Bad.distY[,g],function(v){ifelse(v<2.205e-308,2.205e-308,v)})
      Group.distY[,g] <- Good.distY[,g]+(1-deltaY[g])*Bad.distY[,g]

      Mix.dist[,g] <- prmtrs$X$pi.g[g]*Group.distX[,g]*Group.distY[,g]
    }

    for (g in 1:G){
      zig.hat[,g] <- Mix.dist[,g]/rowSums(Mix.dist)
      vig[,g] <- Good.distX[,g]/Group.distX[,g]
      uig[,g] <- Good.distY[,g]/Group.distY[,g]

    }
    loglikelihood <- sum(log(rowSums(Mix.dist)))

    prmtrs$X$n.g <-prmtrs$Y$n.g <- colSums(zig.hat)
    prmtrs$X$pi.g <- prmtrs$Y$pi.g <- colMeans(zig.hat)
    loglik[i] <- loglikelihood

    deltaX <- colSums(zig.hat*vig)/prmtrs$X$n.g
    deltaY <- colSums(zig.hat*uig)/prmtrs$X$n.g

    a1 <- aitkens(lval = loglik, i = i, eps = tol)
    if (i > 2)
      a2 <- loglik[i - 1]
    for (g in 1:G){
      curr.muX[,g] <- prmtrs$X$mu[,g]
    }

    for (g in 1:G) {
      v1 <- zig.hat[,g]*(vig[,g]*E2.X[,g]+((1-vig[,g])/etaX[g])*E2.tilde.X[,g])
      v2 <- zig.hat[,g]*(vig[,g]*E1.X[,g]+(1-vig[,g])*E1.tilde.X[,g])
      v4 <- zig.hat[,g]*(vig[,g]+((1-vig[,g])/sqrt(etaX[g])))
      A <- sum(v1)
      B <- sum(v2)
      D <- sum(v4)
      mu <- c((B*(v1%*%x)-D*(v4%*%x))/(B*A-D^2))
      cv <- check.val(x,mu,mu.tol)
      if(cv == 0){
        alpha <- c((A*(v4%*%x)-D*(v1%*%x))/(B*A-D^2))
      }else{
        mu <- curr.muX[,g]
        alpha <- c(v4%*%t(t(x)-mu)/B)
      }
      prmtrs$X$mu[,g] <- mu;
      prmtrs$X$alpha[,g] <- alpha;
      prmtrs$X$talpha[g,] <- alpha;

      if (dX==1){
        m1 <- (1/prmtrs$X$n.g[g])*crossprod(sqrt(v1)*(X-prmtrs$X$mu[,g]),sqrt(v1)*(X-prmtrs$X$mu[,g]))
        m2 <- (-1/prmtrs$X$n.g[g])*sum(v4*t(apply(X,1,function(v){(v-prmtrs$X$mu[,g])})))%o%prmtrs$X$alpha[,g]
      }else{
        m1 <- (1/prmtrs$X$n.g[g])*matrix(c(colSums(v1*t(apply(x,1,function(v){(v-prmtrs$X$mu[,g])%o%(v-prmtrs$X$mu[,g])})))),dX,dX)
        m2 <- (-1/prmtrs$X$n.g[g])*colSums(v4*t(apply(x,1,function(v){(v-prmtrs$X$mu[,g])})))%o%prmtrs$X$alpha[,g]
      }
      m3 <- (1/prmtrs$X$n.g[g])*B*prmtrs$X$alpha[,g]%o%prmtrs$X$alpha[,g]


      S <- m1+m2+t(m2)+m3

      prmtrs$X$sig[,,g] <- S
      prmtrs$X$inv.sig[,,g] <- solve(prmtrs$X$sig[,,g])
      prmtrs$X$log.det.sig[g] <- log(det(matrix(prmtrs$X$sig[,,g],dX,dX)))

      #cmstep 2
      v1 <- zig.hat[,g]*(1-vig[,g])
      v2 <- zig.hat[,g]*(1-vig[,g])*E2.tilde.X[,g]
      v3 <- apply(X,1,function(v){(v-prmtrs$X$mu[,g])%*%prmtrs$X$inv.sig[,,g]%*%(v-prmtrs$X$mu[,g])})
      v4 <- apply(X,1,function(v){(v-prmtrs$X$mu[,g])%*%prmtrs$X$inv.sig[,,g]%*%prmtrs$X$alpha[,g]})
      a <- dX*sum(v1)
      b <- c(v1%*%v4)
      c <- c(v2%*%v3)
      sq.eta <- (-b+sqrt(b^2+4*a*c))/(2*a)
      ###############################
      etaX[g] <- max(1,sq.eta^2)
      etaX[g] <- min(etaX[g],etamax)
    }
    ### Y
    curr.muY <- prmtrs$muY
    for (g in 1:G) {
      v1 <- zig.hat[,g]*(uig[,g]*E2.Y[,g]+((1-uig[,g])/etaY[g])*E2.tilde.Y[,g])
      v2 <- zig.hat[,g]*(uig[,g]*E1.Y[,g]+(1-uig[,g])*E1.tilde.Y[,g])
      v4 <- zig.hat[,g]*(uig[,g]+((1-uig[,g])/sqrt(etaY[g])))
      A <- sum(v1)
      B <- sum(v2)
      D <- sum(v4)

      bnum <- crossprod(sqrt(v1)*Xstar,sqrt(v1)*Y)-t(crossprod((v4)%*%(v4%*%Y),Xstar))/B
      bden <- crossprod(sqrt(v1)*Xstar,sqrt(v1)*Xstar)-crossprod(v4%*%(v4%*%Xstar),Xstar)/B
      beta <- solve(bden) %*% (bnum)
      muY <-  Xstar %*% beta
      cv <- check.valY(Y,muY,mu.tol)
      muY[cv==1] <- curr.muY[,,g][cv==1]
      alpha<- (v4%*%(Y-muY))/B
      beta <- solve(t(Xstar) %*% Xstar) %*% t(Xstar) %*% muY
      prmtrs$beta[,,g]<- beta
      prmtrs$muY[,,g] <- muY
      prmtrs$Y$alpha[,g] <- alpha;
      prmtrs$Y$talpha[g,] <- alpha;

      if (dY==1){
        m1 <- (1/prmtrs$X$n.g[g])*crossprod(sqrt(v1)*(Y-prmtrs$muY[,,g]),sqrt(v1)*(Y-prmtrs$muY[,,g]))
        m2 <- (-1/prmtrs$X$n.g[g])*crossprod(sqrt(v4)*(Y-prmtrs$muY[,,g]),sqrt(v4)*prmtrs$Y$alpha[,g])
        m3 <- (1/prmtrs$X$n.g[g])*B*crossprod(prmtrs$Y$alpha[,g],prmtrs$Y$alpha[,g])

        # m1 <- (1/n.g[g])*crossprod(sqrt(v1)*(X-prmtrs$X$mu[,g]),sqrt(v1)*(X-prmtrs$X$mu[,g]))
        # m2 <- (-1/n.g[g])*sum(v4*t(apply(X,1,function(v){(v-prmtrs$X$mu[,g])})))%o%prmtrs$X$alpha[,g]
      }else{
        m1 <- (1/prmtrs$X$n.g[g])*crossprod(sqrt(v1)*(Y-prmtrs$muY[,,g]),sqrt(v1)*(Y-prmtrs$muY[,,g]))
        m2 <- (-1/prmtrs$X$n.g[g])*colSums((v4*(Y-prmtrs$muY[,,g]))%o%prmtrs$Y$alpha[,g])
        m3 <- (1/prmtrs$X$n.g[g]) * B * (prmtrs$Y$alpha[,g] %*% t(prmtrs$Y$alpha[,g]))
      }
      S <- m1+m2+t(m2)+m3
      prmtrs$Y$sig[,,g] <- S
      prmtrs$Y$inv.sig[,,g] <- solve(prmtrs$Y$sig[,,g])
      prmtrs$Y$log.det.sig[g] <- log(det(matrix(prmtrs$Y$sig[,,g],dY,dY)))

      #cmstep 2
      v1 <- zig.hat[,g]*(1-uig[,g])
      v2 <- zig.hat[,g]*(1-uig[,g])*E2.tilde.Y[,g]
      v3 <- diag((Y-prmtrs$muY[,,g])%*%prmtrs$Y$inv.sig[,,g]%*%t(Y-prmtrs$muY[,,g]))
      v4 <- (Y-prmtrs$muY[,,g])%*%prmtrs$Y$inv.sig[,,g]%*%(prmtrs$Y$alpha[,g])

      a <- dY*sum(v1)
      b <- c(v1%*%v4)
      c <- c(v2%*%v3)
      sq.eta <- (-b+sqrt(b^2+4*a*c))/(2*a)
      ###############################
      etaY[g] <- max(1,sq.eta^2)
      etaY[g] <- min(etaY[g],etamax)
    }



    if (i == max.it)
      break
    i <- i + 1
    if (print.iter==T){
      cat(i," ")}
  }
  group  <- apply(zig.hat,1,which.max)
  innerX <- numeric(n)
  innerY <- numeric(n)
  for(k in 1:n){
    innerX[k] <- ifelse(vig[k,group[k]]<0.5,"bad.X","*")
    innerY[k] <- ifelse(uig[k,group[k]]<0.5,"bad.Y","*")
  }

  detection <- data.frame(group=group,detectX=innerX,detectY=innerY)


  AIC   <- -2*loglikelihood + 2*npar
  BIC   <- -2*loglikelihood + npar*log(n)

  z.const    <- (z<10^(-322))*10^(-322)+(z>10^(-322))*z   # vincolo per evitare i NaN nel calcolo di tau*log(tau)
  hard.z     <- (matrix(rep(apply(z,1,max),G),n,G,byrow=F)==z)*1


  return(list(
    Y         = Y,
    X         = X,
    G         = G,
    n         = n,
    salcwm    = salcwm,
    npar      = npar,
    prior     = prmtrs$X$pi.g,
    muX       = prmtrs$X$mu,
    SigmaX    = prmtrs$X$sig,
    alphaX    = prmtrs$X$alpha,
    deltaX    = 1-deltaX,
    etaX      = etaX,
    beta      = prmtrs$beta,
    muY       = prmtrs$muY,
    SigmaY    = prmtrs$Y$sig,
    alphaY    = prmtrs$Y$alpha,
    deltaY    = 1-deltaY,
    etaY      = etaY,
    z         = zig.hat,
    group     = group,
    detection = detection,
    density   = Mix.dist,
    iter.stop = i,
    loglik    = loglikelihood,
    AIC       = AIC,
    BIC       = BIC
  ))

}
