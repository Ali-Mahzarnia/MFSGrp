\name{MFSGrp-package}
\alias{MFSGrp}
\docType{package}
\title{
\packageTitle{MFSGrp}
}
\description{
\packageDescription{MFSGrp}
}
\details{
}
\usage{
MFSGrp(
  Ytrain,
  Xtrain, 
  basisno=5,
  tt, 
  lambda=NULL, 
  alpha=NULL,
  loss="ls",
  part, 
  Xpred=NULL,
  Ypred=NULL, 
  Silence=FALSE, 
  bspline=FALSE, 
  Penalty=NULL, 
  lambdafactor=0.005,
  nfolds=5, 
  predloss="L2",
  eps = 1e-08,
  maxit = 3e+08,
  nlambda=100,
  forcezero=FALSE, 
  forcezeropar=0.001,
  sixplotnum=1, 
  lambdaderivative=NULL,
  nfolder=5,
  nalpha=9,
  nlamder=10,
  lamdermin=1e-9,
  lamdermax=1e-3,
  alphamin=0 ,
  alphamax=1
  )
}
\arguments{
\item{Ytrain}{a vector of scalar response value of size \eqn{n} in the train set.}

\item{Xtrain}{a three dimensional matrix (array) of \eqn{p} functional covariates observed at \eqn{nt} time points and sampled \eqn{n} times in the train set. The size of the array is \eqn{p \times n \times  nt}{p*n*nt}. It is recommended that the time series signals are of the same type. For instance in fMRI data they are all BOLD contrast activities of the human brain that are recorded \eqn{nt} times during a time period \eqn{T} at \eqn{p} Regions of Interests of the Brain Atlas for \eqn{n} human subjects. In another example, in an  economic data set the covariates are time series signals that are all measured in percentage.}

\item{basisno}{the number of bases (of Fourier or B-spline) used in the process of converting the (time series) signals, i.e: observed functional covariates, to functional objects.}

\item{tt}{a vector of  the size \eqn{nt} that includes grid points within 0 to 1. This vector represents the time of observations during the whole period of time \eqn{T} scaled to the interval \eqn{(0,1)}.}

\item{lambda}{ regularization parameter of the (group) lasso penalty term (given alpha=0). By leaving this argument, the program computes the optimized value in a \code{nfolds}-fold cross validation on the train set. The cross validation takes place for values within \code{lambdafactor}*\code{lambdamax} (computed maximum lambda by KKT condition) and \code{lambdamax}, for \code{nlambda} possible values equally spaced in the interval based on the \eqn{log()} of the interval boundaries. if the curve of RMSE vs. log(\code{lambda}) is not U-shaped consider a smaller value for \code{lambdafactor} than the default.}

\item{alpha}{the Elastic Net mixing parameter. For the (group) lasso, \code{alpha} becomes \eqn{0} automatically, for the Ridge penalty it  becomes 1. If it is 1 the penalty becomes ridge automatically. The only way to incorporate a user-defined \code{alpha} value between 0 and \eqn{1} is to set \code{penalty="gelast"}. When \code{penalty="gelast"}, by leaving this argument, the program computes the optimized value in a \code{nfolder}-fold net search for a grid of ten \code{lambda} between the maximum and the minimum \code{lambda} on the train set. The net search cross validation of \code{alpha} takes place for values within \eqn{0} \eqn{<} \code{alphamin} \eqn{<} \code{alphamax} \eqn{<} \eqn{1}, and for \code{nalpha} possible values equally spaced in the interval. The value \code{nfolder} must be greater than 2.}

\item{loss}{inherited from \code{gglasso} package. By default it is set for the least square problem.}
\item{part}{a vector of size \eqn{p} with entries of \code{basisno}. This vector determines the partitions associated to each functional covariate.}
\item{Xpred, Ytrain}{ the test set versions of \code{Xtrain} and \code{Ytrain}.}
\item{Silence}{ by default is \code{"FALSE"}. Given it is \code{"FALSE"}, if in the simulation the population coefficient curves, objects b1,b2,b3,b4,b5 etc. exist they will be displayed by green color along with their estimations. Otherwise, the estimated curves and the number of curves from \eqn{1} to \eqn{p} will be displayed in groups of six figures. There will be \code{sixplotnum} sets of six figures. If \code{sixplotnum = "max"} there will be as many figures as needed to display all of the estimated nonzero curves. At least six curves are displayed if \code{Silence=FALSE}. If in the simulation b1,b2,b3,b4,b5, etc. exist, remove them before running the regression on a real data set. See the example below.}

\item{bspline}{if \code{"FALSE"} the functional conversions and estimations will be processed with Fourier bases. By default is \code{FALSE}.}
 \item{Penalty}{takes four values:
\itemize{ \item{  \code{"glasso"} for (group) lasso.} 
 \item{\code{"ridge"} for ridge.}
 \item{\code{"OLS"} for ordinary least square.}
 \item{\code{"gelast"} for (group) elastic net. See \code{alpha} for more info. }}}
 
\item{loss}{inherited from \code{gglasso} package. By default it is set for the least square \code{"L2"}.}

\item{eps}{inherited from \code{gglasso} package. The threshold of the inner loops. By default it is \code{1e-08}.}

\item{maxit}{inherited from \code{gglasso} package. The maximum iterations of the inner loops.  By default it is \code{3e+08}.}

\item{forcezero}{if \code{"TRUE"}, the estimated coefficient curves the norm of which are less than  \code{forcezeropar} will become zero curves. By default is \code{FALSE}.}

\item{lambdaderivative}{the regularization parameter of the second derivative penalty of the estimated curves. The value \eqn{0} does not penalize the second derivative. By leaving this argument, the program computes the optimized value in a \code{nfolder}-fold net search for a grid of ten \code{lambda} between the maximum and the minimum \code{lambda} on the train set. The net search cross validation of \code{lambdaderivative} takes place for values within \code{lamdermin} \eqn{<} \code{lamdermax}, and for \code{nlamder} possible values equally spaced in the interval based on the \eqn{log()} of the interval boundaries.
If both \code{alpha} and \code{lambdaderivative} need cross validations, first \code{lambdaderivative} is tuned with \code{alpha=0}, then \code{alpha} is regularized with the regularized \code{lambdaderivative}. Both of these processes take place before the regularizations of \code{lambda}. Note that the value \code{nfolder} must be greater than 2.
}
}
\value{
\item{coef}{vector of size \code{basisno}\eqn{*p}. The coefficients of the \eqn{p} estimated functional curves based on \code{basisno} bases.}

\item{predict}{ predicted values associated with \code{Ytest}.}

\item{MSEpredict}{mean squared error of the predictions.}

\item{lambda}{the regularized \code{lambda}.}
}
\examples{

p=35 # number of functional predictor for each observation
n=200 # sample size
nt=500 # number of recorded time points for each functional covariate
# nt will be reduced to 100 after the inner products are computed below

X= array(NaN, c(p, n, nt)); # Brownian motion
for(j in 1:p){
  for (i in 1:n){
    X[j,i,]=cumsum(rnorm(nt,0,1)) }
}

#true nonzero coefs: beta_5, beta_8, beta_11 the rest are zeros
#beta_5(t)=sin(3*pi*t), beta_8(t)=sin(5*pi*t/2) and beta_11(t)=t^2
beta5 = function(t){return(sin(3*pi*t/2))}
beta8 = function(t){return(sin(5*pi*t/2))}
beta11=function(t){return(t^2)}
b5=matrix(0, ncol=1, nrow=nt)
b8=matrix(0, ncol=1, nrow=nt)
b11=matrix(0, ncol=1, nrow=nt)

#evaluate population betas on (0,1) at five hundred time points
for(i in 0:nt){
  j=i
  b5[i]=beta5(j/nt)
  b8[i]=beta8(j/nt)
  b11[i]=beta11(j/nt)
}

#evaluate inner product of Xs and beta 5 and 8 and 11 via Reiman sum
Xb5=matrix(0, ncol=n, nrow=1)
Xb8=matrix(0, ncol=n, nrow=1)
Xb11=matrix(0, ncol=n, nrow=1)

for(j in 1:n){
  Xb5[j]=(X[5,j,] \%*\%b5)/nt
  Xb8[j]=(X[8,j,]\%*\%b8)/nt
  Xb11[j]=(X[11,j,]\%*\%b11)/nt
}
# construct Y
Y=matrix(0, ncol=n, nrow=1)
# standard deviation of the noise term
sd=0.05
# noise term 
eps=matrix(0, ncol=n, nrow=1)
for(n in 1:n){
  eps[, n]=rnorm(1,0,sd)
}
Y=Xb5+Xb8+Xb11+eps
# the algorithm takes care of the intercept in the prediction.
Y=Y+3; #intercept


# make the design matrix (pick every 5 elements) nt becomes 100
X.obs = X[,,(1:100)*nt/100, drop=FALSE]

# observed times scaled to (0,1)
tt=(1:100)/100

# test and train sets (half, half)
trainIndex=sample(1:n, size = round(0.5*n), replace=FALSE)
Ytrain=Y[, trainIndex, drop = FALSE ]
Ytest=Y[, -trainIndex, drop = FALSE ]
Xtrain=X.obs[,trainIndex,, drop=F]
Xtest=X.obs[,-trainIndex,, drop=F]

# The model:
# total 35 functional predictors
# beta_5(t), beta_8(t), beta_11(t) are nonzero, the others are zero

# plot X^1, ..., X^9   (out of total p=35)
par(mfrow=c(3,3))
for(j in 1:9){
  plot(tt,Xtrain[j,1,],type='l', ylim=c(-30,30), main=paste0("j=",j))
  for(k in 2:10)lines(tt, Xtrain[j,k,])
}
par(mfrow=c(1,1))

# load the library 
library(MFSGrp) 

# run
m=17 #basisino
part=rep(m,p) # partition

# green: true beta  (only beta5,beta8,beta11 are the nonzero functions)
# black: estimated beta

# lasso 
# in order to see all figures after the run, use previous plot arrow on Rstudio
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
           Ypred=Ytest, Penalty = "glasso" , bspline=T, sixplotnum="max" , 
           lamdermax=1e-3, lamdermin = 1e-5)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients.
results$lambda # the regularized lambda


# lasso with forcezro 
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
          Ypred=Ytest, Penalty = "glasso" , bspline=T, sixplotnum="max" , 
           lamdermax=1e-3, lamdermin = 1e-5, forcezero=T)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients.


# lasso with forcezro and without penalizing the second derivative
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,
          Ypred=Ytest, Penalty = "glasso" , bspline=T, sixplotnum="max" , 
           lambdaderivative=0, forcezero=T)
sqrt(results$MSEpredict)  # test Root MSE
sum(results$coef==0)/m    # number of zero functional coefficients.


# elastic net with differenrt settings
results=MFSGrp(Ytrain,Xtrain,basisno=m,tt, part=part,Xpred=Xtest,Ypred=Ytest, 
          Penalty = "gelast" , bspline=F, nlamder=15, nalpha=11, nfolder=4 )
sum(results$coef==0)/m
sqrt(results$MSEpredict)


# oracle
active=c(5,8,11)
results=MFSGrp(Y=Ytrain,X=Xtrain[active,,],basisno=m,tt, part=part,Xpred
              =Xtest[active,,],Ypred=Ytest, Penalty ="OLS", bspline = T, 
              sixplotnum = "max" ,lambdaderivative = NULL, lamdermax = 1e-2,
              lamdermin=1e-10, nlamder=20)
sqrt(results$MSEpredict)


#  regression with no shrinkage or sparse penalty, only with second derivative penalty
results=MFSGrp(Ytrain,X=Xtrain,basisno=m,tt, part=part,Xpred=Xtest,Ypred=Ytest, 
          Penalty ="OLS", bspline = T, sixplotnum = 2)
sqrt(results$MSEpredict)
sum(results$coef==0)/m


#  regression with no penalty
results=MFSGrp(Ytrain,X=Xtrain,basisno=m,tt, part=part,Xpred=Xtest,Ypred=Ytest, 
          Penalty ="OLS", bspline = T, sixplotnum = 2, lambdaderivative=0)
sqrt(results$MSEpredict)
sum(results$coef==0)/m


# remove the population curves before data application
remove(b5,b8,b11)
}
\references{
Ali Mahzarnia, Jun Song. "Multivariate functional covariate selection", Submitted in April 2021.

Yang, Y. and Zou, H. (2015), "A Fast Unified Algorithm for
Computing Group-Lasso Penalized Learning Problems," \emph{Statistics and
Computing}. 25(6), 1129-1141.

Friedman, J., Hastie, T., and Tibshirani, R. (2010), "Regularization paths
for generalized linear models via coordinate descent," \emph{Journal of
Statistical Software}, 33, 1.
}
\seealso{
\code{fGMD} package. Install \code{fGMD} in order to have \code{MFSGrp} package work.
}
\author{
Ali Mahzarnia and Jun Song\cr Maintainer: Ali Mahzarnia <amahzarn@uncc.edu>
}
\keyword{Scalar-on-Functional Regression}
\keyword{Elastic net}