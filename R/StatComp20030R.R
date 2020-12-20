#' @title A confidence interval and hypothesis test using R
#' @description A confidence interval and hypothesis test of single normal variance using R
#' @param x  samples
#' @param mu the mean value of the normal distribution
#' @param sigma the variance of the normal distribution
#' @param hypothesis inspection required
#' @param a confidence level
#' @return acceptance or rejection and the confidence interval of estimation
#' @examples
#' \dontrun{
#' A<-rnorm(100,0,2)
#' rnR<-sigma.test(A,mu=0,sigma=2,hypothesis='equal',a=0.05)
#' rnR
#' }
#' @importFrom stats qchisq var
#' @export
sigma2test<-function(x,mu,sigma,hypothesis,a){
  n=length(x)
  s2=var(x)
  
  if(hypothesis=='equal'){
    interl=(n-1)*s2/qchisq(1-(a/2),n-1)
    interr=(n-1)*s2/qchisq(a/2,n-1)
    inter=c(interl,interr)
    c=(n-1)*s2/sigma^2
    if(c>qchisq(a/2,n-1)&&c<qchisq(1-(a/2),n-1))
      AoR='accept'
    else AoR='reject'
    result=list('original hypothesis'=AoR,'interval estimation'=inter)
  }  
  if(hypothesis=='greater'){
    interl=(n-1)*s2/qchisq(1-(a/2),n-1)
    interr=(n-1)*s2/qchisq(a/2,n-1)
    inter=c(interl,interr)
    c=(n-1)*s2/sigma^2
    if(c<qchisq(a,n-1))
      AoR='reject'
    else AoR='accept'
    result=list('original hypothesis'=AoR,'interval estimation'=inter)
  }
  if(hypothesis=='less'){
    interl=0
    interr=(n-1)*s2/qchisq((1-a),n-1)
    inter=c(interl,interr)
    c=(n-1)*s2/sigma^2
    if(c>qchisq(1-a,n-1))
      AoR='reject'
    else AoR='accept'
    result=list('original hypothesis'=AoR,'interval estimation'=inter)
  }
  return(result) 
}

#' @title Random number of Laplacian distribution using R
#' @description Random number of Laplacian distribution using R
#' @param n the number of samples
#' @param beta location parameter of Laplacian distribution
#' @param mu the scale parameter of Laplacian distribution
#' @return a random sample of size \code{n}
#' @examples
#' \dontrun{
#' x<-laplacerandnum(1e4,3,0)
#' hist(x, freq = FALSE,main = "Histogram with the density curve", xlab = "value")
#' }
#' @importFrom stats  runif
#' @export
laplacerandnum<-function(n,beta,mu){
  u<-runif(n)
  x<-numeric(n)
  for (i in 1:n){
    if(u[i]>=0.5)
      x[i]<-mu-beta*log(1-u[i])
    else
      x[i]<-mu+beta*log(u[i])
  }
  return(x)
}



