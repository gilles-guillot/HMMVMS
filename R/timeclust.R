#' @title Model-based clustering for one-dimensional circular data.
#' @description Model-based clustering with two clusters for one-dimensional hour of the day data. Expectation maximization algorithm is implemented that takes into account that the hour of the day data are circular, i.e. 00:00 is the same as 24:00. The circularity is being handled by defining a truncated normal distribution.
#' @param x a vector with the hour of the day data.
#'
#' @return 
#'  \item{loglik }{Final log-likelihood estimate of the EM algorithm.}
#' \item{parameters }{Parameters inferred by the algorithm:
#'     \itemize{
#'       \item pro: Mixing proportion of each distribution.
#'       \item mean: Means of the two clusters.
#'       \item sigma: Sigmas of the two clusters.
#'       \item z: Responsibilities.}}
#'       \item{classification }{Classification of the datapoints to the clusters.}
#' \item{iter }{Number of iterations of the algorithm.}
#' @export
#'
#' @examples
#' x1 <- rnorm(100, 0, 3)
#' x1[x1<0] <- x1[x1<0]+24
#' x2 <- rnorm(80, 9, 2)
#' x <- c(x1,x2)
#' x <- x[sample.int(length(x))]
#' res <- timeclust(x)
#' cat('mixing proportions:',res$parameters$pro,'\nmeans:',
#'     res$parameters$mean,'\nsigmas:',res$parameters$sigma)

timeclust <-
function(x)
{
  quadratic <- function(x, mu, c)
  {
    dist <- x-mu
    if (dist > 12)
      dist <- -24+dist
    else{
      if(dist < -12)
        dist <- 24+dist
    }
    c*dist*dist
  }
  
  cylinderlik <- function(x, mu1, sigma1, mu2, sigma2, tau, int1, int2)
  {
    c1 <- -0.5/sigma1^2; c2 <- -0.5/sigma2^2
    d1 <- tau*(2*pi)^(-1) /sigma1 /int1; d2 <- (1-tau)*(2*pi)^(-1) /sigma2 /int2
    
    suma <- sapply(x, function(x){
      log(d1*exp(quadratic(x,mu1,c1)) + d2*exp(quadratic(x,mu2,c2)))
    })
    return(sum(suma))
  }
  
  responsibilities <- function(x, mu1, sigma1, mu2, sigma2, tau, int1, int2)
  {
    c1 <- -0.5/sigma1^2; c2 <- -0.5/sigma2^2
    d1 <- tau*(2*pi)^(-1) /sigma1 /int1; d2 <- (1-tau)*(2*pi)^(-1) /sigma2 /int2
    
    dat <- sapply(x, function(x) c(d1*exp(quadratic(x,mu1,c1)),
                                   d2*exp(quadratic(x,mu2,c2))) )
    
    mix1 <- apply(dat, 2, function(x) x[1]/(x[1]+x[2]))
    mix2 <- apply(dat, 2, function(x) x[2]/(x[1]+x[2]))
    return(cbind(mix1,mix2,deparse.level=0))
  }
  
  mu_max <- function(x, mu, mixpro, i)
  {
    mu_opt2 <- function(mu, dat, mixpro)
    {
      suma <- apply(cbind(dat, mixpro, deparse.level=0), 1, function(y) {
        dist <- y[1]-mu
        if (dist > 12) {
          dist <- -24+dist
        }else{
          if(dist < -12)
            dist <- 24+dist
        }
        return(y[2]*dist)
      })
      return(sum(suma))
    }
    mu <- nleqslv(mu,mu_opt2,dat=x,mixpro=mixpro[,i])$x
    if (mu>24) {
      mu <- mu-24
    }else {
      if (mu<0)
        mu <- mu+24}
    return(mu)
  }
  
  sigma_max <- function(x, mu, prop, i, diag)
  {
    dat <- cbind(x,prop[,i],deparse.level=0)
    sigmas <- apply(dat, 1, function(x) {
      dist <- x[1]-mu
      if (dist > 12)
        dist <- -24+dist
      else{
        if(dist < -12)
          dist <- 24+dist
      }
      return(dist^2*x[2])})
    sigmasq <- sum(sigmas)/sum(prop[,i])
    return(sqrt(sigmasq))
  }
  
  # Initialization of the parameters to be estimated.
  tau <- 0.5
  mu1 <- mean(head(x,100))
  mu2 <- max(head(x,100))
  sigma1 <- 1
  sigma2 <- 1
  
  int1 <- 1-2*pnorm(mu1-12,mu1,sigma1)
  int2 <- 1-2*pnorm(mu2-12,mu2,sigma2)
  
  cllk <- rep(NA, 1000)
  cllk[1] <- 0
  cllk[2] <- cylinderlik(x, mu1, sigma1, mu2, sigma2, tau, int1, int2)
  k<-2
  
  # loop
  while(abs(cllk[k]-cllk[k-1]) >= 0.00001) {
    
    # E step  
    mixpro <- responsibilities(x, mu1, sigma1, mu2, sigma2, tau, int1, int2)
    
    # M step  
    tau <- sum(mixpro[,1])/length(mixpro[,1])
    
    mu1 <- mu_max(x, mu1, mixpro, 1)
    mu2 <- mu_max(x, mu2, mixpro, 2)
    
    sigma1 <- sigma_max(x, mu1, mixpro, 1)
    sigma2 <- sigma_max(x, mu2, mixpro, 2)
    
    int1 <- 1-2*pnorm(mu1-12,mu1,sigma1)
    int2 <- 1-2*pnorm(mu2-12,mu2,sigma2)
    
    cllk[k+1] <- cylinderlik(x, mu1, sigma1, mu2, sigma2, tau, int1, int2)
    k <- k+1
  }
  
  return(list(loglik=cllk[k], parameters=list(pro=c(tau, 1-tau), mean=c(mu1,mu2),
                                              sigma=c(sigma1,sigma2), z=mixpro),
              classification=apply(mixpro, 1, function(x) ifelse(x[1]>x[2],1,2)), iter=(k-1)))
}
