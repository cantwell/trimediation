#' CounterfactuaL Disparity Measure (CDM) of the inverse probability-
#' weighted (IPW) marginal structural model.
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
#' @return : (vector) a length-2 list with the CDM with RR-interaction and
#'   RD-interaction.
#' @export
#' @seealso \code{\link{cde_f}}
#' @seealso \code{\link{stran_f}}
#' @seealso \code{\link{g_est_f}}
#' @seealso \code{\link{tml_f}}
#' @references Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M.
#'   (2016). Mediation Analysis for Health Disparities Research. American
#'   Journal of Epidemiology, 184(4), 315–324.
#'   https://doi.org/10.1093/aje/kwv329
#' @author Ashley Naimi \email{ashley.naimi@pitt.edu}
#' @author Cantwell Carson \email{carsonc@gmail.com}
#' @examples see cde_example.R
ipw_f<-function(outcome,mY=NULL,famY=NULL,
                mediator,mM,famM,
                exposure,mX,famX,
                data,interactionRR,interactionRD){

  X_propensity<-glm(mX,data=data,family=famX)$fitted.values
  M_propensity<-glm(mM,data=data,family=famM)$fitted.values

  den<-(X_propensity*exposure + (1-X_propensity)*(1-exposure))*(M_propensity*mediator + (1-M_propensity)*(1-mediator))
  num<-(mean(exposure)*exposure + (1-mean(exposure))*(1-exposure))*(mean(mediator)*mediator + (1-mean(mediator))*(1-mediator))

  sw=num/den

  if(interactionRR==T){
    RR<-exp(coef(glm(outcome~exposure+mediator+exposure:mediator,data=data,family=poisson(link="log"),weights=sw))["exposure"])
    names(RR)<-"risk ratio"
  } else{
    RR<-exp(coef(glm(outcome~exposure+mediator,data=data,family=poisson(link="log"),weights=sw))["exposure"])
    names(RR)<-"risk ratio"
  }

  if(interactionRD==T){
    RD<-coef(glm(outcome~exposure+mediator+exposure:mediator,data=data,family=gaussian(link="identity"),weights=sw))["exposure"]
    names(RD)<-"risk difference"
  } else{
    RD<-coef(glm(outcome~exposure+mediator,data=data,family=gaussian(link="identity"),weights=sw))["exposure"]
    names(RD)<-"risk difference"
  }

  ipw<-cbind(RD,RR);row.names(ipw)<-"IPW"

  return(ipw)

}
