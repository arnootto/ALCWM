# The SALCWM package
An R package to implement the methodology described in *Mixture of multivariate linear asymmetric Laplace regressions with multiple asymmetric Laplace covariates* (2025).

To install the package, use the following code in R
```{r}
#install.packages("devtools")
library(devtools)
install_github("arnootto/SALCWM")
```
## Example
Code to reproduce Example 7.1: Australian Institute of Sport data
```{r}
library(sn)
data("ais")
data=ais
subset=data[,c("sex","RCC","WCC","BMI","SSF","Bfat","LBM")]
est1=SALCWM(Y=cbind(data$RCC,data$WCC),X=cbind(data$BMI,data$SSF,data$Bfat,data$LBM),G=2,tol=1e-5, max.it=2000,initialization = "mclust")
est2=CSALCWM(Y=cbind(data$RCC,data$WCC),X=cbind(data$BMI,data$SSF,data$Bfat,data$LBM),G=2, tol=1e-5,max.it=2000,initialization = "mclust")
```
