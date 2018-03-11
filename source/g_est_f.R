#' CounterfactuaL Disparity Measure (CDM) of the G-estimation model.
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
#' @seealso \code{\link{ipw_f}}
#' @seealso \code{\link{stran_f}}
#' @seealso \code{\link{tml_f}}
#' @references Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M.
#'   (2016). Mediation Analysis for Health Disparities Research. American
#'   Journal of Epidemiology, 184(4), 315–324.
#'   https://doi.org/10.1093/aje/kwv329
#' @author Ashley Naimi \email{ashley.naimi@pitt.edu}
#' @author Cantwell Carson \email{carsonc@gmail.com}
#' @examples see cde_example.R
g_est_f<-function(outcome,mY=NULL,famY=NULL,
              mediator,mM,famM,
              exposure,mX,famX,
              data,interactionRR=T,interactionRD=T){

  X_propensity<-glm(mX,data=data,family=famX)$fitted.values
  M_propensity<-glm(mM,data=data,family=famM)$fitted.values

  if(interactionRD==T){
    mod1<-matrix(coef(lm(y~x+I(m-M_propensity)+x:I(m-M_propensity)+c_my,data=data))[c("I(m - M_propensity)","x:I(m - M_propensity)")])
    y_tilde<-outcome - matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1
    RD<-coef(lm(y_tilde~I(x-X_propensity)+c_xy,data=data))["I(x - X_propensity)"]
  } else{
    mod1<-matrix(coef(lm(y~x+I(m-M_propensity)+c_my,data=data))[c("I(m - M_propensity)")])
    y_tilde<-outcome - matrix(mediator,ncol=1)%*%mod1
    RD<-coef(lm(y_tilde~I(x-X_propensity)+c_xy,data=data))["I(x - X_propensity)"]
  }

  if(interactionRR==T){
    mod1<-matrix(coef(glm(y~x+I(m-M_propensity)+x:I(m-M_propensity)+c_xy+c_my,family=poisson,data=data))[c("I(m - M_propensity)","x:I(m - M_propensity)")])
    y_tilde<-outcome*exp( - matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1)
    RR<-exp(coef(glm(y_tilde~I(x-X_propensity)+c_xy,family=quasipoisson,data=data))["I(x - X_propensity)"])
  } else{
    mod1<-matrix(coef(glm(y~x+I(m-M_propensity)+c_xy+c_my,family=poisson,data=data))[c("I(m - M_propensity)")])
    y_tilde<-outcome*exp( - matrix(mediator,ncol=1)%*%mod1)
    RR<-exp(coef(lm(y_tilde~I(x-X_propensity)+c_xy,family=quasipoisson,data=data))["I(x - X_propensity)"])
  }

  names(RD)<-"Risk Difference"
  names(RR)<-"Risk Ratio"
  g_est<-cbind(RD,RR);row.names(g_est)<-"G Estim"

  return(g_est)
}
