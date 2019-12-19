#' @title Model-based clustering for two-dimensional linear-circular data.
#' @description Model-based clustering with two clusters for two-dimensional data, where one of the dimensions is the hour of the day. Expectation maximization algorithm is implemented that takes into account that the hour of the day data are circular, i.e. 00:00 is the same as 24:00. The circularity is being handled by defining a truncated normal distribution.
#'
#' @param x a matrix with the linear data in the first column and the circular data in the second one.
#' @param diag a boolean indicating if the covariance matrix of the joint distribution is diagonal, i.e. the linear and the circular variables are independent.
#'
#' @return 
#' loglik Final log-likelihood estimate of the EM algorithm.
#' parameters Parameters inferred by the algorithm:
#' pro: Mixing proportion of each distribution.
#'  mean: Means of the two clusters.
#'      sigma: Covariance matrices of the two clusters.
#'       z: Responsibilities.
#'       classification: Classification of the datapoints to the clusters.
#'       iter: Number of iterations of the algorithm.
#' @export
#'
#' @examples 
#' require(mixtools)
#' x1 <- mixtools::rmvnorm(n= 250, mu= c(3.5,12), sigma= matrix(c(1,0,0,4),nrow=2,ncol=2))
#' x2 <- mixtools::rmvnorm(n= 350, mu= c(3,22), sigma= matrix(c(1,0,0,6),nrow=2,ncol=2))
#' x2[x2[,2]>24, 2] <- x2[x2[,2]>24, 2] - 24
#' x <- rbind(x1,x2,deparse.level=0)
#' x <- x[sample.int(nrow(x)),]
#' res <- circlust(x)
#' cat('mixing proportions:',res$parameters$pro,'\nmean of cluster 1:\n',res$parameters$mean[,1],
#'     '\ncovariance matrix of cluster 1:\n',res$parameters$sigma[,,1],
#'     '\nmean of cluster 2:\n',res$parameters$mean[,2],
#'     '\ncovariance matrix of cluster 2:\n',res$parameters$sigma[,,2])
#' 
#' 
circlust <-
function(x, diag=F)
{
  quadratic <- function(x, mu, a, b, c)
  {
    dist1 <- x[1]-mu[1]
    dist2 <- x[2]-mu[2]
    if (dist2 > 12)
      dist2 <- -24+dist2
    else{
      if(dist2 < -12)
        dist2 <- 24+dist2
    }
    a*dist1*dist1 + b*dist1*dist2 + c*dist2*dist2
  }
  
  cylinderlik <- function(x, mu1, sigma1, mu2, sigma2, tau, int1, int2)
  {
    det1 <- sigma1[1,1]*sigma1[2,2]-sigma1[2,1]^2
    det2 <- sigma2[1,1]*sigma2[2,2]-sigma2[2,1]^2
    a1 <- -0.5*sigma1[2,2]/det1; a2 <- -0.5*sigma2[2,2]/det2
    b1 <- sigma1[2,1]/det1; b2 <- sigma2[2,1]/det2
    c1 <- -0.5*sigma1[1,1]/det1; c2 <- -0.5*sigma2[1,1]/det2
    d1 <- tau*(2*pi)^(-1)*det1^(-0.5) /int1; d2 <- (1-tau)*(2*pi)^(-1)*det2^(-0.5) /int2
    
    suma <- apply(x, 1, function(x){
      log(d1*exp(quadratic(x,mu1,a1,b1,c1)) + d2*exp(quadratic(x,mu2,a2,b2,c2)))
    })
    return(sum(suma))
  }
  
  responsibilities <- function(x, mu1, sigma1, mu2, sigma2, tau, int1, int2)
  {
    det1 <- sigma1[1,1]*sigma1[2,2]-sigma1[2,1]^2
    det2 <- sigma2[1,1]*sigma2[2,2]-sigma2[2,1]^2
    a1 <- -0.5*sigma1[2,2]/det1; a2 <- -0.5*sigma2[2,2]/det2
    b1 <- sigma1[2,1]/det1; b2 <- sigma2[2,1]/det2
    c1 <- -0.5*sigma1[1,1]/det1; c2 <- -0.5*sigma2[1,1]/det2
    d1 <- tau*(2*pi)^(-1)*det1^(-0.5) /int1; d2 <- (1-tau)*(2*pi)^(-1)*det2^(-0.5) /int2
    
    dat <- apply(x, 1, function(x) c(d1*exp(quadratic(x,mu1,a1,b1,c1)),
                                     d2*exp(quadratic(x,mu2,a2,b2,c2))) )
    
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
    mu <- c(sum(x[,1] * mixpro[,i])/sum(mixpro[,i]), nleqslv(mu[2],mu_opt2,dat=x[,2],mixpro=mixpro[,i])$x)
    if (mu[2]>24) {
      mu[2] <- mu[2]-24
    }else {
      if (mu[2]<0)
        mu[2] <- mu[2]+24}
    return(mu)
  }
  
  sigma_max <- function(x, mu, prop, i, diag)
  {
    Tau <- sum(prop[,i])
    dat <- cbind(x,prop[,i],deparse.level=0)
    sigmas <- apply(dat, 1, function(x) {
      dist2 <- x[2]-mu[2]
      if (dist2 > 12)
        dist2 <- -24+dist2
      else{
        if(dist2 < -12)
          dist2 <- 24+dist2
      }
      c((x[1]-mu[1])^2*x[3] , (x[1]-mu[1])*dist2*x[3] , dist2*dist2*x[3])})
    
    sigmas <- rowSums(sigmas)
    if (diag)
      sigma <- matrix(c(sigmas[1]/Tau, 0, 0, sigmas[3]/Tau), nrow=2, ncol=2)
    else
      sigma <- matrix(c(sigmas[1]/Tau, sigmas[2]/Tau, sigmas[2]/Tau, sigmas[3]/Tau), nrow=2, ncol=2)
    return(sigma)
  }
  
  # Initialization of the parameters to be estimated.
  tau <- 0.5
  mu1 <- c(mean(head(x[,1],100)),mean(head(x[,2],100)))
  mu2 <- c(mean(head(x[,1],100)),max(head(x[,2],100)))
  sigma1 <- matrix(c(1,0,0,1), nrow=2, ncol=2)
  sigma2 <- matrix(c(1,0,0,1), nrow=2, ncol=2)
  if (diag){
    int1 <- 1-2*pnorm(mu1[2]-12,mu1[2],sigma1[2,2])
    int2 <- 1-2*pnorm(mu2[2]-12,mu2[2],sigma2[2,2])
  }else{
    int1 <- 1-2*pmvnorm(upper=c(Inf,mu1[2]-12),mean=mu1,sigma=sigma1)
    int2 <- 1-2*pmvnorm(upper=c(Inf,mu2[2]-12),mean=mu2,sigma=sigma2)}
  
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
    
    sigma1 <- sigma_max(x, mu1, mixpro, 1, diag)
    sigma2 <- sigma_max(x, mu2, mixpro, 2, diag)
    
    if (diag){
      int1 <- 1-2*pnorm(mu1[2]-12,mu1[2],sigma1[2,2])
      int2 <- 1-2*pnorm(mu2[2]-12,mu2[2],sigma2[2,2])
    }else{
      int1 <- 1-2*pmvnorm(upper=c(Inf,mu1[2]-12),mean=mu1,sigma=sigma1)
      int2 <- 1-2*pmvnorm(upper=c(Inf,mu2[2]-12),mean=mu2,sigma=sigma2)}
    
    cllk[k+1] <- cylinderlik(x, mu1, sigma1, mu2, sigma2, tau, int1, int2)
    k <- k+1
  }
  
  return(list(loglik=cllk[k], parameters=list(pro=c(tau, 1-tau), mean=cbind(mu1,mu2),
                                              sigma=array(c(sigma1,sigma2),dim=c(2,2,2)), z=mixpro),
              classification=apply(mixpro, 1, function(x) ifelse(x[1]>x[2],1,2)), iter=(k-1)))
}
