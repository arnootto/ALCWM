#' Mixture of multivariate linear asymmetric Laplace regressions with multiple asymmetric Laplace covariates
#'
#' This package fits the shifted asymmetric Laplace cluster-weighted model and the contaminated shifted asymmetric Laplace cluster-weighted model via an EM-algorithm, as described in "Mixture of multivariate linear asymmetric Laplace regressions with multiple asymmetric Laplace covariates" (2025).
#'
#'@section Key Functions:
#' \itemize{
#'   \item \code{\link{SALCWM}}: Fits the Shifted Asymmetric Laplace Cluster-Weighted Model.
#'   \item \code{\link{CSALCWM}}: Fits the Contaminated Shifted Asymmetric Laplace Cluster-Weighted Model.
#'   \item \code{\link{dSAL}}: Evaluates the density of the Shifted Asymmetric Laplace distribution.
#'   \item \code{\link{dCSAL}}: Evaluates the density of the Contaminated Shifted Asymmetric Laplace distribution.
#' }
#'
#'@examples
#'# Fit the SALCWM and CSALCWM models with G = 2 components to the AIS dataset
#'library(ALCWM)
#'library(sn)
#'data("ais")
#'est1=SALCWM(Y=cbind(ais$RCC,ais$WCC),X=cbind(ais$BMI,ais$SSF,ais$Bfat,ais$LBM),G=2,tol=1e-5, max.it=2000,initialization = "mclust")
#'est2=CSALCWM(Y=cbind(ais$RCC,ais$WCC),X=cbind(ais$BMI,ais$SSF,ais$Bfat,ais$LBM),G=2, tol=1e-5,max.it=2000,initialization = "mclust")
#'
#' @docType package
#' @name ALCWM
#' @aliases ALCWM-package
#' @author Arno Otto <arno.otto@up.ac.za>, Andriëtte  Bekker, Antonio Punzo, Johan Ferreira, Cristina Tortora.
#' @keywords package
NULL
