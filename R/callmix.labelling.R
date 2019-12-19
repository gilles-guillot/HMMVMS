#' Classification of positions as fishing or steaming, using speed data.
#'
#' @param x  a vector of the speed data
#' @param nclust the number of clusters
#' @param algo 'depmix' or 'mclust'
#' @details By using the methods provided in the \code{algo} argument, model-based clustering is performed, with the specified number of clusters. The clusters are then assigned to fishing of steaming using the following labelling scheme: Fishing corresponds to the combination of components that has a joint empirical variance smaller than that of the low-speed states estimated in a 2-component model. If there is more than one combination that reduces the variance, we choose the one with a joint empirical mean closest to that of the low-speed states estimates in a 2-component model.
#' @return A list of 3 elements:
#' `states` A vector of the classified fishing or steaming states.
#' `loglik` Final log-likelihood estimate of the EM algorithm.
#' `parameters` A list of parameters dependent on the model.
#'
#' @examples data("DanishTrips")
#' res <- callmix.labelling(DanishTrips$speed,3)
#' @export
callmix.labelling <-
function(x,nclust,algo = 'depmix')
{
  ##
  call.mclust <- function(x,nclust)
  {
    ## independent mixture of two Gaussians
    try.mclust = try(Mclust(data=x,G=nclust),silent=TRUE)
    return(try.mclust)
  }
  ##
  call.depmix <- function(x,nclust)
  {
    x = data.frame(x=x)
    model = depmix(response = x ~ 1, data=x, nstates=nclust)
    try.depmix = try(fit(model,verbose=FALSE),silent=TRUE)
    return(try.depmix)
  }
  
  
  ###
  callmix2 <- function(x,nclust,algo)
  {
    ## 0 are handled before the call to this function
    isa = !is.na(x)
    
    classification <- rep(NA,length(x))
    mean <- rep(NA,nclust)
    sd <- rep(NA,nclust)
    posterior <- matrix(NA, nrow = length(x), ncol = nclust)
    logLik <- NA
    parameters <- NA
    
    ## If there isn't enough pings, we don't try to infer
    if (sum(isa) > 10)
    { 
      if (algo == "mclust")
      {
        y <- x[isa]
        res <- call.mclust(x = y, nclust = nclust)
        if(class(res) == 'try-error' | is.null(res)) {return(NULL)}
        classification[isa] <- res$classification
        
        posterior[isa,] <- res$z
        
        mean <- res$parameters$mean
        if(res$parameters$variance$modelName=='V'){
          sd <-sqrt(res$parameters$variance$sigmasq) }
        if(res$parameters$variance$modelName=='E'){
          sd <-rep(sqrt(res$parameters$variance$sigmasq),
                   nclust) }
        parameters <- res$parameters
        logLik <- res$loglik
      }
      if (algo == "depmix")
      {
        y <- x
        res <- call.depmix(x = y, nclust = nclust)
        if (class(res) == 'try-error' | is.null(res)) {return(NULL)}
        classification <- res@posterior[,1]
        posterior <- data.matrix(res@posterior[,2:(nclust+1)])
        for(iclust in 1:nclust)
        {
          mean[iclust] <- res@response[[iclust]][[1]]@parameters$coefficients
          sd[iclust] <- res@response[[iclust]][[1]]@parameters$sd
        }
        params <- getpars(res)
        init <- params[1:nclust]
        transition <- matrix(params[seq(nclust+1,nclust*(nclust+1))],ncol = nclust, nrow = nclust)
        parameters = list(init = init, transition = transition, mean = mean, sd = sd)
        logLik <- logLik(res)
      }
      ## Label switching
      oorder <- order(mean)
      mean <- mean[oorder]
      sd <- sd[oorder]
      oorder.inv = order(oorder)
      classification <- oorder.inv[classification]
      posterior <- posterior[,oorder]
      ## Label switching on parameters
      if (algo == "depmix") { 
        parameters$init <- parameters$init[oorder] 
        names(parameters$init) <- names(parameters$init)[oorder]
        parameters$transition <- parameters$transition[oorder,oorder]
        parameters$mean <- parameters$mean[oorder]
        parameters$sd <- parameters$sd[oorder]
      }
      if (algo == "mclust") {
        parameters$pro <- parameters$pro[oorder]
        parameters$mean <- parameters$mean[oorder]
        names(parameters$mean) <- names(parameters$mean)[oorder]
        if(parameters$variance$modelName=='V'){
          parameters$variance$sigmasq <- parameters$variance$sigmasq[oorder]
          parameters$variance$scale   <- parameters$variance$scale[oorder]
        }
      }
      res <- list(classification = classification, mean = mean, sd = sd, 
                  posterior = posterior, logLik = logLik, parameters = parameters)
    }
    return(res)
  }
  
  if (nclust == 2) ## Meaning 2 states fishing/steaming
  {
    res.callmix <- callmix2(x = x,nclust = nclust, algo = algo)
    if (is.null(res.callmix)) {return(NULL)}
    classification = res.callmix$classification == 1
    logLik = res.callmix$logLik
    parameters = res.callmix$parameters
    mean <- res.callmix$mean
    sd <- res.callmix$sd
  }
  if (nclust > 2) ## Meaning n states with labelling
  {
    res.to.label <- callmix2(x = x, nclust = nclust, algo = algo)
    res.ref      <- callmix2(x = x, nclust = 2     , algo = algo)
    if (is.null(res.to.label) || is.null(res.ref) || all(is.na(res.to.label$classification))) {return(NULL)}
    
    ## Initialization
    classification = rep(NA,length(x))
    
    logLik = res.to.label$logLik
    parameters = res.to.label$parameters
    
    classification.to.label <- res.to.label$classification
    mean.1clust <- res.to.label$mean
    sd.1clust <- res.to.label$sd
    posterior <- res.to.label$posterior
    
    mean.fish <- res.ref$mean[1]
    sd.fish <- res.ref$sd[1]
    
    class.no.na <- classification.to.label[!is.na(classification.to.label)]
    clusts <- sort(unique(class.no.na))
    ## Computing all permutations
    comb.clust <- expand.grid( as.data.frame( matrix( rep(0:1,length(clusts)), nrow=2  ) ) ) == 1
    ## Removing no fishing combination
    comb.clust <- comb.clust[apply(comb.clust,1,any),,drop = FALSE]
    ## Compute the mixing proportions of the resulting mixture distribution and mean/sd of combinations
    mix.prop <- matrix(0,nrow = dim(comb.clust)[1], ncol = dim(comb.clust)[2])
    mean.comb <- rep(NA,dim(comb.clust)[1])
    sd.comb   <- rep(NA,dim(comb.clust)[1])
    ind.clust <- matrix(FALSE, nrow = dim(comb.clust)[1], ncol = length(class.no.na))
    for (icomb in seq_len(dim(comb.clust)[1]))
    {
      ind.clust[icomb,] <- is.element(class.no.na,clusts[comb.clust[icomb,]])
      if (sum(comb.clust[icomb,])>1) {
        post.prob.comb <- posterior[ind.clust[icomb,],comb.clust[icomb,], drop = FALSE]
        post.prob.comb <- post.prob.comb / apply(post.prob.comb,1,sum)
        mix.prop[icomb,comb.clust[icomb,]] <- apply(post.prob.comb,2,sum) / sum(ind.clust[icomb,])
      } else {
        mix.prop[icomb,comb.clust[icomb,]] <- 1
      }
      mean.comb[icomb] <- sum(mix.prop[icomb,] * mean.1clust)
      sd.comb[icomb]   <- sqrt( sum( mix.prop[icomb,] * ( (mean.1clust - mean.comb[icomb])^2 +
                                                            sd.1clust^2) ) )
    }
    if (any(sd.comb < sd.fish)) {
      ind.best <- which(sd.comb<sd.fish)[which.min(abs(mean.comb[sd.comb<sd.fish] - mean.fish))]
    } else {
      ind.best <- which.min(abs(mean.comb - mean.fish))
    }
    classification[!is.na(classification.to.label)] <- ind.clust[ind.best,]
    ## Computing non-fishing mixing proportions
    post.prob.nonF <- posterior[!ind.clust[ind.best,],!comb.clust[ind.best,], drop = FALSE]
    post.prob.nonF <- post.prob.nonF / apply(post.prob.nonF,1,sum)
    mix.prop.nonF  <- apply(post.prob.nonF,2,sum) / sum(ind.clust[ind.best,])
    mean.nonF <- sum(mix.prop.nonF * mean.1clust[!comb.clust[ind.best,]])
    sd.nonF <- sqrt( sum( mix.prop.nonF * ( (mean.1clust[!comb.clust[ind.best,]] - mean.nonF)^2 +
                                              sd.1clust[!comb.clust[ind.best,]]^2) ) )
    mean <- c(fishing=mean.comb[ind.best], steaming=mean.nonF)
    sd <- c(fishing=sd.comb[ind.best], steaming=sd.nonF)
    parameters <- c(meanFS=mean, sdFS=sd, parameters)
  }
  ## Transform TRUE/FALSE in Fishing/Steaming
  ind.class <- classification
  classification[ ind.class] <- 'F'
  classification[!ind.class] <- 'S'
  res <- list(states = classification, loglik = logLik, parameters = parameters)
  return(res)
}
