#' @title Classification of VMS positions as fishing or steaming
#' @description Classifies positions as fishing or steaming using different techniques with VMS data. Possible techniques are: model-based clustering on speed data (with packages mclust or mixtools), model-based clustering with dependences on speed data (with the depmix package), model-based clustering on speed and hour of the day data (with the custom functions timeclust and circlust that into account that the hour of the day is circular data).
#' @details The labeling of the clusters required to perform the classification is performed by the function(s) \code{\link{callmix.labelling}}, \code{classify.mixtools}, \code{\link{classify.circular}}.
#' @param data a data.frame with the VMS data. First column should contain time in POSIXlt format. Second column should contain the speed data. Columns 3 and 4 should contain the vessel id and the trip id respectively. 
#' @param methods a list of vectors with first element the method name and second element the number of clusters (if necessary). Different methods can be: depmix, mclust, mixtools, timeclust, circlust. If NULL it calculates the estimates from depmix and mclust with 2 and 3 clusters.
#' @param analysis a character with the type of analysis, 'vessel' for vessel by vessel analysis, 'trip' for trip by trip analysis or 'all' to perform the analysis on the whole data set.
#' @param runs an integer denoting the number of different runs of each method. (The inference doesn't guarantee a global optimum, so the maximum likelihood solution is chosen among different runs.)
#' 
#' @return A list of 2 elements:
#' a data.frame with the input data and the resulting classification for each of the methods provided by the methods argument of the function.
#' a list of lists containing the parameters infered by the relevant models. (Each list element, for each trip/ vessel of the analysis contains a list of parameters for each of the methods) In the depmix case, it is the outcome of the \code{getpars} command. In the mclust case, it is the parameters list of the output.
#' @export
#'
#' @examples 
#' res <- classifyTrips(DanishTrips, analysis='all')
#' misclassification <- function(v, s)
#' {
#'   d <- 0
#'   for (i in 1:length(v))
#'     if (v[i] != s[i])
#'       d <- d+1
#'     return(d/length(v)*100)
#' }
#' misclassification(res$states[[6]],DanishTrips$state_seq)

classifyTrips <-
  function(data, methods=list(c('depmix',3)), analysis='vessel', runs=3)
{
  
  classify.normalmixEM <- function(fitted)
  {
    fi <- which.min(fitted$sigma)
    fa <- which.max(fitted$mu)
    if (fi==fa)
      fi <- which.min(fitted$sigma[-fi])
    classification <- apply(fitted$posterior, 1, which.max)
    state_seq.est <- ifelse(classification==fi, 'F', 'S')
    return(state_seq.est)
  }
  
  estimate.sequence <- function(frame, method, nclust=NULL)
  {
    if (method=='depmix'){
      result <- callmix.labelling(frame[,'speed'], nclust, method)
    }else if (method=='mclust'){
      result <- callmix.labelling(frame[,'speed'], nclust, method)
    }else if (method=='mixtools'){
      fitted <- try(normalmixEM(frame[,'speed'],k=nclust), silent=T)
      if (class(fitted) != 'try-error'){
        result <- list(states=classify.normalmixEM(fitted), loglik=fitted$loglik, parameters=
                         list(lambda=fitted$lambda, mu=fitted$mu, sigma=fitted$sigma))
      }
    }else if (method=='timeclust'){
      fitted <- try(classify.circular(frame[,'speed'], frame[,'time'], method),silent=T)
      if (class(fitted) != 'try-error'){
        result <- list(states=fitted$classification, loglik=fitted$loglik, parameters=fitted$parameters)
      }
    }else if (method=='circlust'){
      fitted <- try(classify.circular(frame[,'speed'], frame[,'time'], method),silent=T)
      if (class(fitted) != 'try-error'){
        result <- list(states=fitted$classification, loglik=fitted$loglik, parameters=fitted$parameters)
      }
    }else{
      print('nonexistent method')
    }
    if (!exists('result')){
      result <- NULL
    }
    return(result)
  }
  
  estimate.sequence.multiple.runs <- function(small, method, nclust, runs){
    statesmat <- matrix(nrow=nrow(small), ncol=runs)
    llk <- numeric()
    param <- list()
    for (run in 1:runs){
      temp <- estimate.sequence(small, method, nclust)
      if (!is.null(temp)){
        statesmat[,run] <-  temp[['states']]
        llk[run] <- temp[['loglik']]
        param[[run]] <- temp[['parameters']]
      }
    }
    if (sum(is.na(statesmat[1,])) < runs/2){
      res <- list(states=statesmat[, which.max(llk)], parameters=param[[which.max(llk)]])
    }else{
      res <- NULL
    }
    return(res)
  }
  
  small.estimate <- function(frame, methods, runs)
  {
    small <- data.frame(matrix(nrow=nrow(frame),ncol=length(methods)),stringsAsFactors=F)
    parameters <- rep(list(NULL),length(methods))
    for (i in 1:length(methods)){
      temp <- estimate.sequence.multiple.runs(frame,methods[[i]][1],as.integer(methods[[i]][2]), runs)
      if (!is.null(temp)){
        small[,i] <- temp[['states']]
        parameters[[i]] <- temp[['parameters']]
      }
    }
    names(parameters) <- lapply(methods, function(x) paste(x[1],x[2],sep=''))
    return(list(frame=small, parameters=parameters))
  }
  
  
  if (sum(is.na(data)) != 0){
    print('Data should be cleaned. i.e. remove NAs.')
    return(NULL)
  }
  dates <- as.vector(as.character(data[,1]))
  if ('POSIXlt' %in% class(dates)){
    print('The time column must be in the POSIXlt format.')
    return(NULL)}
  dates <- as.timeDate(dates)
  time <- hour(dates)+minute(dates)/60
  speed <- as.vector(as.numeric(data[,2]))
  if (class(speed)!='numeric'){
    print('The speed column must be in numeric format.')
    return(NULL)}
  vessels_vec <- as.vector(as.character(data[,3]))
  if (class(vessels_vec)!='character'){
    print('The vessels column must be in character format.')
    return(NULL)}
  trips_vec <- as.vector(as.character(data[,4]))
  if (class(trips_vec)!='character'){
    print('The trips column must be in character format.')
    return(NULL)}
  
  frame <- data.frame(time=time, speed=speed, vessel=vessels_vec, trip=trips_vec, stringsAsFactors=F)
  colnames(frame) <- c('time','speed','vessel','trip')
  indexes <- frame[,'speed']>1e-4
  frame <- frame[indexes,]
  vessels <- unique(frame[,'vessel'])
  trips <- unique(frame[,'trip'])
  
  extended.frame <- data.frame(matrix(nrow=nrow(frame),ncol=length(methods)), stringsAsFactors=F)
  colnames(extended.frame) <- lapply(methods, function(x) paste(x[1],x[2],sep=''))
  
  if (analysis=='all'){
    temp <- small.estimate(frame, methods, runs)
    parameters <- temp[['parameters']]
    extended.frame[,] <- temp[['frame']]
    
  }else if (analysis=='trip'){
    parameters <- list()
    for (i in 1:length(trips)){
      cat('Classifying trip',trips[i],'\n')
      cond <- frame[,'trip']==trips[i]
      temp <- small.estimate(frame[cond,], methods, runs)
      parameters[[i]] <- temp[['parameters']]
      extended.frame[cond, ] <- temp[['frame']]
    }
    names(parameters) <- trips
  }else{
    parameters <- list()
    for (i in 1:length(vessels)){
      cat('Classifying vessel',vessels[i],'\n')
      cond <- frame[,'vessel']==vessels[i]
      temp <- small.estimate(frame[cond,], methods, runs)
      parameters[[i]] <- temp[['parameters']]
      extended.frame[cond, ] <- temp[['frame']]
    }
    names(parameters) <- vessels
  }
  temp <- data.frame(matrix('S', nrow=nrow(data), ncol=length(methods)), stringsAsFactors=F)
  colnames(temp) <- colnames(extended.frame)
  temp[indexes,] <- extended.frame
  frame.res <- data.frame(data, temp, stringsAsFactors=F)
  
  return(list(states=frame.res, parameters=parameters))
}
