#' Simulated datasets for demonstration
#' 
#' The dataset was generated based on the proposed model \eqn{h(t; \boldsymbol{Z}_i, {X}_i(\cdot))=h_{0}(t-t_{i,j-1}) \exp \left(\eta_{ij}\right)},
#' where \eqn{h_0(\cdot)} is the baseline hazard function generated from a Weibull distribution. \eqn{\eta_{ij} = \bm{\alpha}^{\top}\boldsymbol{Z}_i +\int_{t_{i, j-1}}^{t}{X}_{i}(s)\beta(s)ds + v_{ij}}. 
#' \eqn{\bm{\alpha}} is the fixed effect parameter associated with the time-invariant covariates \eqn{\boldsymbol{Z}_i}, 
#' and \eqn{\beta(t)} is a time-varying coefficient that captures the effect of functional predictor \eqn{X_{i}(t)} on the hazard rate of recurrent events.
#' The simulated dataset is organized into two data frames:
#' a survival data frame (\code{sdat}) and a functional data frame (\code{fdat}).
#' The variables in each data frame are listed below:
#' 
#'
#' @usage data(simDat)
#' @source Simulated data
#' @docType data
#' @name simDat
#' @aliases sdat fdat simData
#' @format A list with two data frame: 
#' \describe{
#'   \item{sdat}{Survival data; a data frame with xxx rows and xxx variables:}
#'   \describe{
#'     \item{id}{Subjects identification}
#'     \item{event}{A sequence of the number of events per subject}
#'     \item{t_start}{Event starting time}
#'     \item{t_end}{Event end time}
#'     \item{censoring_time}{Event censoring time}
#'     \item{status}{Event status: \code{status=1} if event is observed and \code{status=0} if event is censored}
#'     \item{z1}{A univariate scalar covariates. Can be extended to multiple scalar covariates}
#'   }
#'   \item{fdat}{Functional data; a data frame with xxx rows and xxx variables:}
#'   \describe{
#'     \item{id}{Subjects identification}
#'     \item{time}{Time points for each longitudinal measurement}
#'     \item{x}{Longitudinal measurements at distinct time points}
#'   }
#'  }
NULL
