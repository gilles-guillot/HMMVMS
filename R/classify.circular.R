#' Classification of positions as fishing or steaming, using speed and hour of the day data.
#' @description Classification of positions as fishing or steaming, using speed and hour of the day data, with model-based clustering.
#' @param speed  a vector of speed data.
#' @param time a vector of hour of the day data.
#' @param method 'timeclust' or 'circlust' 
#' @param diag  a boolean to be passed to  'circlust' 
#' @details The classification is performed in two stages. First, 
#' a clustering with two clusters is performed on the speed data. 
#' The cluster with the high speed mean is classified as steaming. 
#' The rest of the positions are used as an input for the `timeclust``
#' or `circlust` methods. From the resulting classification, 
#' the positions that have a time mean closer to 9 are classified as fishing and the rest as steaming.
#' @return A list of 3 elements: 
#' `classification` A vector of the classified fishing or steaming states.
#' `loglik` Final log-likelihood estimate of the EM algorithm. 
#' `parameters` A list of parameters dependent on the model.
#' @export
#'
#' @examples dates <- as.timeDate(DanishTrips$DATES)
#' time <- hour(dates)+minute(dates)/60
#' res <- classify.circular(DanishTrips$speed, time, 'timeclust')
#' misclassification <- function(v, s)
#' {
#'   d <- 0
#'   for (i in 1:length(v))
#'     if (v[i] != s[i])
#'       d <- d+1
#'     return(d/length(v)*100)
#' }
#' misclassification(res$classification, DanishTrips$state_seq)
classify.circular <-
function(speed, time, method, diag=F)
{
  state.est <- vector()
  modelo=try(Mclust(speed,G=2,modelNames='V'),silent=T)
  if (class(modelo) != 'try-error'){
    goodstate=which(modelo$parameters$mean == min(modelo$parameters$mean))
    indexes=which(modelo$classification == goodstate)
  }else{
    modelo=try(normalmixEM(speed),silent=T)
    if (class(modelo) != 'try-error'){
      goodstate=which(modelo$mu == min(modelo$mu))
      classification=apply(modelo$posterior, 1, function(x) ifelse(x[1]>x[2],1,2))
      indexes=which(classification == goodstate)
    }
  }
  if (class(modelo) != 'try-error'){
    newtime=time[indexes]
    if (method=='circlust'){
      newdata=speed[indexes]
      x=cbind(newdata,newtime,deparse.level=0)  
      modelo2=try(circlust(x,diag),silent=TRUE)
      if (class(modelo2) != 'try-error'){
        if (abs(modelo2$parameters$mean[2,1]-9) < abs(modelo2$parameters$mean[2,2]-9))
          check <- T
      }
    }else{
      modelo2=try(timeclust(newtime),silent=TRUE)
      if (class(modelo2) != 'try-error'){
        if (abs(modelo2$parameters$mean[1]-9) < abs(modelo2$parameters$mean[2]-9))
          check <- T
      }
    }
    if (class(modelo2) != 'try-error'){
      if (check)
        indexes2=which(modelo2$classification==1)
      else
        indexes2=which(modelo2$classification==2)
      state.temp=rep('S',length(indexes))
      state.temp[indexes2]='F'
      state.est=rep('S',length(speed))
      state.est[indexes]=state.temp
      loglik <- modelo2$loglik
    }
  }
  return(list(classification=state.est, loglik=loglik, parameters=list(model1=modelo$parameters,
                                                                       model2=modelo2$parameters[1:3])))
}
