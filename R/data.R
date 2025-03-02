#' Simulated dataset for demonstration
#'
#' The simulated dataset is organized into two data frames:
#' a survival data frame and a functional data frame.
#' The variables in each data frame are listed below:
#' 
#' @description Survival data \code{sdat}
#' \describe{
#'   \item{id}{subjects identification}
#'   \item{event}{sequence of the number of events per subject}
#'   \item{t_start}{event starting time}
#'   \item{t_end}{event end time}
#'   \item{censoring_time}{event censoring time}
#'   \item{status}{event status: status=1 if event observed and status=0 if event censored}
#'   \item{z1}{a univariate scalar covariates. Can be extended to multiple scalar covariates}
#' }
#' @description Functional data \code{fdat}
#' \describe{
#'   \item{id}{subjects identification}
#'   \item{time}{time points for each longitudinal measurement}
#'   \item{x}{longitudinal measurements at distinct time points}
#' }
#'
#' @usage data(simDat)
#' @docType data
#' @name simDat
#' @rdname simDat
#' @format Two data frames ...
NULL
