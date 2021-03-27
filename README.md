# MFSGrp
This package runs the Group Elastic Net (including lasso, ridge, elastic net, and ordinary least square) regression with scalar response values, and observed functional covariates. In addition, it penalizes the curvature of the output by implementing a penalty on the second derivative of the estimated coefficient curves. A part of this package uses the [fGMD](https://github.com/Ali-Mahzarnia/fGMD) package that is created exclusively for this package. The [fGMD](https://github.com/Ali-Mahzarnia/fGMD) package is a heavily modified version of the [gglasso](https://github.com/cran/gglasso). The features added to the original package: mixing parameter alpha and its net search cross validation, curvature penalizing for functional regression and its net search cross-validations, optimized Fortran core function to speed up the curvature penalization updates, progress reports and time estimations. In order to have this package work first install [fGMD](https://github.com/Ali-Mahzarnia/fGMD). The [fGMD](https://github.com/Ali-Mahzarnia/fGMD) package does not work independently from this package and it does not interfere with the functions of the [gglasso](https://github.com/cran/gglasso) package due to slight name differences.

## Installation:
  ## First install [fGMD](https://github.com/Ali-Mahzarnia/fGMD):
You can install `fGMD` from [GitHub](https://github.com/Ali-Mahzarnia/fGMD) with the R code:
```R
install.packages("https://github.com/Ali-Mahzarnia/fGMD/raw/master/fGMD_1.0.tar.gz",  repos = NULL, type="source")
```

## Then install MFSGrp:
You can install `MFSGrp` from [GitHub](https://github.com/Ali-Mahzarnia/MFSGrp) with the R code:
```R
install.packages("https://github.com/Ali-Mahzarnia/MFSGrp/raw/main/MFSGrp_1.0.tar.gz",  repos = NULL, type="source")
```

## Manual and examples:
After installations you can pull up the manual that includes a simulations example
```R
??MFSGrp
```
## Main refrence
Ali Mahzarnia, Jun Song. "Multivariate functional covariate selection", Submitted in March 2021.
