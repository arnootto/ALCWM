# The ALCWM package
An R package to implement the methodology described in *Mixtures of multivariate linear asymmetric Laplace regressions with multiple asymmetric Laplace covariates* (2025).

To install the package, use the following code in R
```{r}
#install.packages("devtools")
library(devtools)
install_github("arnootto/ALCWM")
```
## Example
Code to reproduce Example 7.1: Australian Institute of Sport data
```{r}
# Fit the SALCWM and CSALCWM models with G = 2 components to the AIS dataset
library(ALCWM)
library(sn)
data("ais")
est1=ml.SALCWM(Y=cbind(ais$RCC,ais$WCC),X=cbind(ais$BMI,ais$SSF,ais$Bfat,ais$LBM),G=2,tol=1e-5, max.it=2000,initialization = "mclust")
est2=ml.cSALCWM(Y=cbind(ais$RCC,ais$WCC),X=cbind(ais$BMI,ais$SSF,ais$Bfat,ais$LBM),G=2, tol=1e-5,max.it=2000,initialization = "mclust")
```
