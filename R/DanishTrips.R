#' @title  Selected variables from the VMS data from 10 Danish trips.
#' @description Time, speed, vessel id, trip id and "true" state sequence from the VMS data collected from 10 Danish trips.
#' @format  A data frame with 671 observations on the following 5 variables.
#' \describe{
#'   \item{\code{DATES}}{a POSIXct}
#'   \item{\code{speed}}{a numeric vector}
#'   \item{\code{vessels_vec}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6}}
#'   \item{\code{trips_vec}}{a factor with levels \code{1} \code{10} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9}}
#'   \item{\code{state_seq}}{a factor with levels \code{F} \code{S}}
#' }
#' @source Bastardie et al. (2010) Detailed mapping of fishing effort and landings by coupling fishing logbooks with satellite-recorded vessel geo-location. Fisheries Research, 106(1) :41-53, 2010.
#'  