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

#' CounterfactuaL Disparity Measure (CDM) of the Targeted Minimum Loss-based
#' Estimation (TMLE) model.
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
#' @references Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M.
#'   (2016). Mediation Analysis for Health Disparities Research. American
#'   Journal of Epidemiology, 184(4), 315–324.
#'   https://doi.org/10.1093/aje/kwv329
#' @author Ashley Naimi \email{ashley.naimi@pitt.edu}
#' @author Cantwell Carson \email{carsonc@gmail.com}
#' @examples see cde_example.R
tml_f<-function(outcome,mY=NULL,famY=NULL,
                mediator,mM=NULL,famM=NULL,
                exposure,mX=NULL,famX=NULL,
                data,interactionRR=T,interactionRD=T){
  mY<-glm(y~c_xy+c_my,data=data,subset=x==1&m==0,family=binomial(link="logit"))
  mM<-glm(m~c_my,data=data,subset=x==1,family=binomial(link="logit"))
  mX<-glm(x~c_xy,data=data,family=binomial(link="logit"))
  newdata1<-data;newdata1$x<-1;newdata1$m<-0;
  pY<-predict(mY,newdata=newdata1,type="response")
  pM<-predict(mM,newdata=newdata1,type="response")
  pX<-predict(mX,newdata=newdata1,type="response")
  data$cc2 <- as.numeric(data$x==1&data$m==0)/((1-pM)*(pX))
  q2_star<-glm(y~-1+cc2,data=data,offset=log(pY/(1-pY)),family=binomial(link="logit"))$fitted.values
  mQ1<-glm(q2_star ~ c_xy,data=data,subset=x==1,family=quasibinomial(link="logit"))
  q1<-predict(mQ1,newdata=newdata1,type="response")
  data$cc1 <- as.numeric(data$x==1)/pX
  mu10<-glm(q2_star ~ -1 + cc1, data=data,family=quasibinomial(link="logit"),offset=log(q1/(1-q1)))$fitted.values

  mY<-glm(y~c_xy+c_my,data=data,subset=x==0&m==0,family=binomial(link="logit"))
  mM<-glm(m~c_my,data=data,subset=x==0,family=binomial(link="logit"))
  mX<-glm(x~c_xy,data=data,family=binomial(link="logit"))
  newdata1<-data;newdata1$x<-0;newdata1$m<-0;
  pY<-predict(mY,newdata=newdata1,type="response")
  pM<-predict(mM,newdata=newdata1,type="response")
  pX<-predict(mX,newdata=newdata1,type="response")
  data$cc2 <- as.numeric(data$x==0&data$m==0)/((1-pM)*(pX))
  q2_star<-glm(y~-1+cc2,data=data,offset=log(pY/(1-pY)),family=binomial(link="logit"))$fitted.values
  mQ1<-glm(q2_star ~ c_xy,data=data,subset=x==0,family=quasibinomial(link="logit"))
  q1<-predict(mQ1,newdata=newdata1,type="response")
  data$cc1 <- as.numeric(data$x==0)/pX
  mu00<-glm(q2_star ~ -1 + cc1, data=data,family=quasibinomial(link="logit"),offset=log(q1/(1-q1)))$fitted.values

  RD<-mean(mu10-mu00)
  RR<-mean(mu10)/mean(mu00)
  names(RD)<-"Risk Difference"
  names(RR)<-"Risk Ratio"
  tml<-cbind(RD,RR);row.names(tml)<-"TMLE"

  return(tml)

}

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

  args = list(outcome,mY,famY, mediator,mM,famM, exposure,mX,famX, data,interactionRR=T,interactionRD=T)

  estimators = list(ipw_f, stran_f, g_est_f, tml_f)

  cde_vals <- data.frame()

  for (estimator in estimators){
    cde_vals <- rbind(cde_vals, do.call(estimator, args))
  }

  return(cde_vals)

}
