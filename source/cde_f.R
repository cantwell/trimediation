#' CounterfactuaL Disparity Measure (CDM) model comparison
#'
#' @param outcome : (datatype) Description of structure. What it represents.
#' @param mY : (datatype) Description of structure. What it represents.
#' @param famY : (datatype) Description of structure. What it represents.
#' @param mediator : (datatype) Description of structure. What it represents.
#' @param mM : (datatype) Description of structure. What it represents.
#' @param famM : (datatype) Description of structure. What it represents.
#' @param exposure : (datatype) Description of structure. What it represents.
#' @param mX : (datatype) Description of structure. What it represents.
#' @param famX : (datatype) Description of structure. What it represents.
#' @param data : (datatype) Description of structure. What it represents.
#' @param interactionRR : (datatype) Description of structure. What it
#'   represents.
#' @param interactionRD : (datatype) Description of structure. What it
#'   represents.
#'
#' @return : (data frame) a (2, 4) data frame with the RR-interaction and
#'   RD-interaction for each of the models included in the analysis.
#' @export
#' @seealso \code{\link{ipw_f.ipw_f}}
#' @seealso \code{\link{stran_f.stran_f}}
#' @seealso \code{\link{g_est_f.g_est_f}}
#' @seealso \code{\link{tml_f.tml_f}}
#' @references Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M.
#'   (2016). Mediation Analysis for Health Disparities Research. American
#'   Journal of Epidemiology, 184(4), 315–324.
#'   https://doi.org/10.1093/aje/kwv329
#' @author Ashley Naimi \email{ashley.naimi@pitt.edu}
#' @author Cantwell Carson \email{carsonc@gmail.com}
#' @examples see cde_example.R
cde_f<-function(outcome,mY=NULL,famY=NULL,
              mediator,mM,famM,
              exposure,mX,famX,
              data,interactionRR=T,interactionRD=T){

  source("R/ipw_f.R")
  source("R/stran_f.R")
  source("R/g_est_f.R")
  source("R/tml_f.R")

  args = list(outcome,mY,famY, mediator,mM,famM, exposure,mX,famX, data,interactionRR=T,interactionRD=T)

  estimators = list(ipw_f, stran_f, g_est_f, tml_f)

  cde_vals <- data.frame()

  for (estimator in estimators){
    cde_vals <- rbind(cde_vals, do.call(estimator, args))
  }

  return(cde_vals)

}
