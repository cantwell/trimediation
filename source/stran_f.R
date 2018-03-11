#' CounterfactuaL Disparity Measure (CDM) of the structural
#' transformation model.
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
#' @seealso \code{\link{g_est_f}}
#' @seealso \code{\link{tml_f}}
#' @references Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M.
#'   (2016). Mediation Analysis for Health Disparities Research. American
#'   Journal of Epidemiology, 184(4), 315–324.
#'   https://doi.org/10.1093/aje/kwv329
#' @author Ashley Naimi \email{ashley.naimi@pitt.edu}
#' @author Cantwell Carson \email{carsonc@gmail.com}
#' @examples see cde_example.R
stran_f<-function(outcome,mY=NULL,famY=NULL,
              mediator,mM=NULL,famM=NULL,
              exposure,mX=NULL,famX=NULL,
              data,interactionRR=T,interactionRD=T){

  if(interactionRD==T){
    mod1<-matrix(coef(lm(y~x+m+x:m+c_xy+c_my,data=data))[c("m","x:m")])
    y_tilde<-outcome - matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1
    RD<-coef(lm(y_tilde~exposure+c_xy,data=data))["exposure"]
  } else{
    mod1<-matrix(coef(lm(y~x+m+c_xy+c_my,data=data))[c("m")])
    y_tilde<-outcome - matrix(mediator,ncol=1)%*%mod1
    RD<-coef(lm(y_tilde~exposure+c_xy,data=data))["exposure"]
  }

  if(interactionRR==T){
    mod1<-matrix(coef(glm(y~x+m+x:m+c_xy+c_my,data=data,family=poisson))[c("m","x:m")])
    y_tilde<-outcome*exp(- matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1)
    RR<-exp(coef(glm(y_tilde~exposure+c_xy,data=data,family=quasipoisson))["exposure"])
  } else{
    mod1<-matrix(coef(glm(y~x+m+c_xy+c_my,data=data,family=poisson))[c("m")])
    y_tilde<-outcome*exp( - matrix(mediator,ncol=1)%*%mod1)
    RR<-exp(coef(glm(y_tilde~exposure+c_xy,data=data,family=poisson))["exposure"])
  }

  names(RD)<-"Risk Difference"
  names(RR)<-"Risk Ratio"
  stran<-cbind(RD,RR);row.names(stran)<-"Struct Transf"

  return(stran)

}
