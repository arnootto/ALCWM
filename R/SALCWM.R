#' shifted asymmetric Laplace cluster-weighted model.
#'
#' Maximum likelihood fitting of the shifted asymmetric Laplace (SAL) cluster-weighted model (CWM) using an EM algorithm.
#'
#' @param Y A matrix or vector containing the response variables.
#' @param X A matrix or vector containing the predictor variables.
#' @param G Number of clusters (components) to fit. Default is 2.
#' @param tol Convergence tolerance for the EM algorithm. Default is 0.001.
#' @param max.it Maximum number of iterations for the EM algorithm. Default is 1000.
#' @param initialization Method for initializing the posterior probabilities. Options are
#'   \code{"mclust"} (default), \code{"kmeans"}, \code{"random.soft"}, \code{"random.hard"},
#'   or \code{"manual"}.
#' @param print.iter Logical. If \code{TRUE}, prints the iteration number during fitting.
#'   Default is \code{TRUE}.
#' @param seed Optional seed for random initialization methods. Default is \code{NULL}.
#' @param start.z Optional matrix of initial posterior probabilities for
#'   \code{initialization = "manual"}. Default is \code{NULL}.
#'
#'' @return A list containing the following components:
#'   \itemize{
#'     \item \code{Y}: The input response matrix.
#'     \item \code{X}: The input predictor matrix.
#'     \item \code{G}: Number of clusters.
#'     \item \code{n}: Number of observations.
#'     \item \code{npar}: Number of model parameters.
#'     \item \code{prior}: Estimated prior probabilities for each cluster.
#'     \item \code{prmtrs}: List of estimated parameters, including \code{X} and \code{Y}
#'       components (means, covariances, skewness parameters).
#'     \item \code{muX}: Estimated mean vectors for \code{X} in each cluster.
#'     \item \code{SigmaX}: Estimated covariance matrices for \code{X} in each cluster.
#'     \item \code{alphaX}: Estimated skewness parameters for \code{X} in each cluster.
#'     \item \code{beta}: Estimated regression coefficients for \code{Y} given \code{X} in
#'       each cluster.
#'     \item \code{muY}: Estimated mean vectors for \code{Y} in each cluster.
#'     \item \code{SigmaY}: Estimated covariance matrices for \code{Y} in each cluster.
#'     \item \code{alphaY}: Estimated skewness parameters for \code{Y} in each cluster.
#'     \item \code{z}: Posterior probabilities for each observation in each cluster.
#'     \item \code{group}: Cluster assignments based on maximum posterior probability.
#'     \item \code{density}: Densities for each observation in each cluster.
#'     \item \code{iter.stop}: Number of iterations until convergence.
#'     \item \code{loglik}: Final log-likelihood value.
#'     \item \code{AIC}: Akaike Information Criterion.
#'     \item \code{BIC}: Bayesian Information Criterion.
#'   }
#' @import mclust
#'
#' @export
SALCWM <- function(Y,
                   X,
                   G=2,
                   tol=0.001,
                   max.it=1000,
                   initialization = "mclust",
                   print.iter=T,
                   seed=NULL,
                   start.z=NULL
){
  if(is.vector(Y))
    Y <- matrix(Y,ncol=1)
  if(is.vector(X))
    X <- matrix(X,ncol=1)
  if(any(is.na(X)) | any(is.na(Y)))
    stop('No NAs allowed.')

   dY <- ncol(Y)
   dX <- ncol(X)
   x <- X
   n  <- nrow(X)

   z  <- array(0,c(n,G),dimnames=list(1:n,paste("comp.",1:G,sep="")))

   npar <- G-1 + G * (dX + dX * (dX + 1)/2 + dX) + G * (dY * (dX + 1) + dY * (dY + 1)/2 + dY)

   prior  <- numeric(G)

   # X Parameters definition
   muX    <- array(0,c(dX,G),dimnames=list(paste("X.",1:dX,sep=""),paste("comp.",1:G,sep="")))
   SigmaX <- array(0,c(dX,dX,G),dimnames=list(paste("X.",1:dX,sep=""),paste("X.",1:dX,sep=""),paste("comp.",1:G,sep="")))
   alphaX <- array(0,c(dX,G),dimnames=list(paste("X.",1:dX,sep=""),paste("comp.",1:G,sep="")))

   # Y Parameters definition
   beta    <- array(0,c(dX+1,dY,G),dimnames=list(c("Intercept",paste("Slope X.",1:dX,sep="")),paste("Y.",1:dY,sep=""),paste("comp.",1:G,sep="")))
   muY     <- array(0,c(n,dY,G),dimnames=list(1:n,paste("Y.",1:dY,sep=""),paste("comp.",1:G,sep="")))
   SigmaY  <- array(0,c(dY,dY,G),dimnames=list(paste("Y.",1:dY,sep=""),paste("Y.",1:dY,sep=""),paste("comp.",1:G,sep="")))
   alphaY     <- array(0,c(dY,G),dimnames=list(paste("Y.",1:dY,sep=""),paste("comp.",1:G,sep="")))

   Xstar <- cbind(rep(1,n),X)

   # Distribution
   dens  <- array(0,c(n,G),dimnames=list(1:n,paste("comp.",1:G,sep="")))
   densX  <- array(0,c(n,G),dimnames=list(1:n,paste("comp.",1:G,sep="")))
   densY  <- array(0,c(n,G),dimnames=list(1:n,paste("comp.",1:G,sep="")))

   E1.X <- E2.X <- E1.Y <- E2.Y <- matrix(NA,n,G)
   prmtrs <-NULL

   #posterior initialization
   if(G>1){

     if(initialization=="mclust"){
       D <- cbind(Y,X)                 # the complete set of data
       model <- Mclust(D,G=G,modelNames="VVV", verbose = print.iter)
       z <- model$z
       cluster <-model$classification
     }


     if(initialization=="kmeans"){
       D         <- cbind(Y,X)                 # the complete set of data
       model  <- kmeans(x=D,centers=G)      # clusters on D
       cluster <- model$cluster
       z         <- clas(model$cluster,G)
     }

     if(initialization=="random.soft"){
       if(!is.null(seed)) set.seed(seed)
       z  <- array(runif(n*G),c(n,G)) # soft posterior probabilities (no-normalized) (n x k)
       z  <- z/rowSums(z)             # soft posterior probabilities (n x k)
       cluster  <- apply(z,1,which.max)


     }

     if(initialization=="random.hard"){
       if(!is.null(seed)) set.seed(seed)
       z  <- t(rmultinom(n, size = 1, prob=rep(1/G,G)))  # hard posterior probabilities (n x k)
       cluster <- apply(z,1,which.max)

     }

     if(initialization=="manual"){ # z.start can be both soft and hard initialization
       z  <- start.z      # posterior probabilities (n x k) no-normalized
       cluster  <- apply(z,1,which.max)

     }

   } else{z <- matrix(rep(1,n),nrow=n,ncol=1)
   cluster  <- matrix(rep(1,n),nrow=n,ncol=1)
   }

   #initialize parameters based on intialized posteriors
   init.prmtrsX<-prmtrs.init(x=X,cluster=cluster,G=G,start=1)
   init.prmtrsY<-prmtrs.init(x=Y,cluster=cluster,G=G,start=1)
   prmtrs$X<-init.prmtrsX
   prmtrs$Y<-init.prmtrsY
   prmtrs$beta <- beta
   prmtrs$muY <- array(rep(t(prmtrs$Y$mu), each = n), dim = c(n, nrow(prmtrs$Y$mu), ncol(prmtrs$Y$mu)))

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
     prior <- colMeans(z)
     n.g <- colSums(z)
     for (g in 1:G){
       e1e2 <- e1e2f(p=dX,x=X,mu=prmtrs$X$mu[,g],alpha=prmtrs$X$alpha[,g],inv.sig=prmtrs$X$inv.sig[,,g])
       E1.X[,g] <- e1e2$E1
       E2.X[,g] <- e1e2$E2
     }
     curr.muX <- prmtrs$X$mu
     for (g in 1:G) {
       v1 <- z[,g]*(E2.X[,g])
       v2 <- z[,g]*(E1.X[,g])
       v4 <- z[,g]
       A <- sum(v1)
       B <- sum(v2)
       D <- sum(v4)
       mu <- c((B*(v1%*%X)-D*(v4%*%X))/(B*A-D^2))
       cv <- check.val(X,mu)
       if(cv == 0){
         alpha <- c((A*(v4%*%X)-D*(v1%*%X))/(B*A-D^2))
       }else{
         mu <- curr.muX[,g]
         alpha <- c(v4%*%t(t(x)-mu)/B)
       }
       prmtrs$X$mu[,g] <- mu;
       prmtrs$X$alpha[,g] <- alpha;
       prmtrs$X$talpha[g,] <- alpha;
       if (dX==1){
         m1 <- (1/n.g[g])*crossprod(sqrt(v1)*(X-prmtrs$X$mu[,g]),sqrt(v1)*(X-prmtrs$X$mu[,g]))
         m2 <- (-1/n.g[g])*sum(v4*t(apply(X,1,function(v){(v-prmtrs$X$mu[,g])})))%o%prmtrs$X$alpha[,g]
       }else{
         m1 <- (1/n.g[g])*matrix(c(colSums(v1*t(apply(x,1,function(v){(v-prmtrs$X$mu[,g])%o%(v-prmtrs$X$mu[,g])})))),dX,dX)
         m2 <- (-1/n.g[g])*colSums(v4*t(apply(x,1,function(v){(v-prmtrs$X$mu[,g])})))%o%prmtrs$X$alpha[,g]
       }
       m3 <- (1/n.g[g])*B*prmtrs$X$alpha[,g]%o%prmtrs$X$alpha[,g]


       S <- m1+m2+t(m2)+m3
       prmtrs$X$sig[,,g] <- S
       prmtrs$X$inv.sig[,,g] <- solve(prmtrs$X$sig[,,g], tol=1e-20)
       prmtrs$X$log.det.sig[g] <- log(det(matrix(prmtrs$X$sig[,,g],dX,dX)))
     }


     ####for Y
     for (g in 1:G){
       e1e2 <- e1e2fY(p=dY,x=Y,mu=prmtrs$muY[,,g],alpha=prmtrs$Y$alpha[,g],inv.sig=prmtrs$Y$inv.sig[,,g])
       E1.Y[,g] <- e1e2$E1
       E2.Y[,g] <- e1e2$E2
     }
     curr.muY <- prmtrs$muY
     for (g in 1:G) {
       v1 <- z[,g]*(E2.Y[,g])
       v2 <- z[,g]*(E1.Y[,g])
       v4 <- z[,g]
       A <- sum(v1)
       B <- sum(v2)
       D <- sum(v4)
       bnum <- crossprod(sqrt(v1)*Xstar,sqrt(v1)*Y)-t(crossprod((v4)%*%(v4%*%Y),Xstar))/B
       bden <- crossprod(sqrt(v1)*Xstar,sqrt(v1)*Xstar)-crossprod(v4%*%(v4%*%Xstar),Xstar)/B
       beta <- solve(bden, tol=1e-20) %*% (bnum)
       muY <-  Xstar %*% beta
       cv <- check.valY(Y,muY)
       muY[cv==1] <- curr.muY[,,g][cv==1]
       alpha<- (v4%*%(Y-muY))/B
       beta <- solve(t(Xstar) %*% Xstar) %*% t(Xstar) %*% muY
       prmtrs$beta[,,g]<- beta
       prmtrs$muY[,,g] <- muY
       prmtrs$Y$alpha[,g] <- alpha;
       prmtrs$Y$talpha[g,] <- alpha;
       if (dY==1){
         m1 <- (1/n.g[g])*crossprod(sqrt(v1)*(Y-prmtrs$muY[,,g]),sqrt(v1)*(Y-prmtrs$muY[,,g]))
         m2 <- (-1/n.g[g])*crossprod(sqrt(v4)*(Y-prmtrs$muY[,,g]),sqrt(v4)*prmtrs$Y$alpha[,g])
         m3 <- (1/n.g[g])*B*crossprod(prmtrs$Y$alpha[,g],prmtrs$Y$alpha[,g])
       }else{
         m1 <- (1/n.g[g])*crossprod(sqrt(v1)*(Y-prmtrs$muY[,,g]),sqrt(v1)*(Y-prmtrs$muY[,,g]))
         m2 <- (-1/n.g[g])*colSums((v4*(Y-prmtrs$muY[,,g]))%o%prmtrs$Y$alpha[,g])
         m3 <- (1/n.g[g]) * B * (prmtrs$Y$alpha[,g] %*% t(prmtrs$Y$alpha[,g]))
       }
       S <- m1+m2+t(m2)+m3
       prmtrs$Y$sig[,,g] <- S
       prmtrs$Y$inv.sig[,,g] <- solve(prmtrs$Y$sig[,,g], tol=1e-20)
       prmtrs$Y$log.det.sig[g] <- log(det(matrix(prmtrs$Y$sig[,,g],dY,dY)))
     }
     for(g in 1:G){
       dens[,g] <- prior[g]*dSAL(X,mu=prmtrs$X$mu[,g],alpha=prmtrs$X$alpha[,g],sigma=prmtrs$X$sig[,,g])*dSALY(Y,mu = prmtrs$muY[,,g], alpha = prmtrs$Y$alpha[,g], sigma= prmtrs$Y$sig[,,g])
       densX[,g] <- prior[g]*dSAL(X,mu=prmtrs$X$mu[,g],alpha=prmtrs$X$alpha[,g],sigma=prmtrs$X$sig[,,g])
       densY[,g] <- prior[g]*dSALY(Y,mu = prmtrs$muY[,,g], alpha = prmtrs$Y$alpha[,g], sigma= prmtrs$Y$sig[,,g])
     }
     loglikelihood <- sum(log(rowSums(dens)))
     loglik[i] <- loglikelihood
     loglikX[i] <- sum(log(rowSums(densX)))
     loglikY[i] <- sum(log(rowSums(densY)))
     is.unsorted(loglikX)
     is.unsorted(loglikY)

     for(g in 1:G){
       z[,g] <- dens[,g]/rowSums(dens)
     }
     if (print.iter==T){
       cat(i," ")}
     a1 <- aitkens(lval = loglik, i = i, eps = tol)
     if (i > 2)
       a2 <- loglik[i - 1]

     if (i == max.it)
       break
     i <- i + 1


   }
   group  <- apply(z,1,which.max)
   AIC   <- -2*loglikelihood + 2*npar
   BIC   <- -2*loglikelihood + npar*log(n)

   z.const    <- (z<10^(-322))*10^(-322)+(z>10^(-322))*z   # vincolo per evitare i NaN nel calcolo di tau*log(tau)
   hard.z     <- (matrix(rep(apply(z,1,max),G),n,G,byrow=F)==z)*1
   return(list(
     Y         = Y,
     X         = X,
     G         = G,
     n         = n,
     npar      = npar,
     prior     = prior,
     prmtrs    = prmtrs,
     muX       = prmtrs$X$mu,
     SigmaX    = prmtrs$X$sig,
     alphaX    = prmtrs$X$alpha,
     beta      = prmtrs$beta,
     muY       = prmtrs$muY,
     SigmaY    = prmtrs$Y$sig,
     alphaY    = prmtrs$Y$alpha,
     z         = z,
     group     = group,
     density   = dens,
     iter.stop = i,
     loglik    = loglikelihood,
     AIC       = AIC,
     BIC       = BIC
   ))

}

