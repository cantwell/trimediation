## 1) revise last three functions using cde_ipw as template -- DONE
## 2) push everything to GH
## 3) implement bootstrap for CI
## 4) test functons in applied datasets

# library(zeallot)

## INTERNAL FUNCTIONS
#  get_propensity
#' This is an an internal function called \code{get_propensity()}
#'
#' This returns the propensities for exposure and mediator
#'
#' @keywords internal
#'
#' @param data : (datatype) Description of structure. What it represents.
#' @param mX : (datatype) Description of structure. What it represents.
#' @param mM : (datatype) Description of structure. What it represents.
#' @param famX : (datatype) Description of structure. What it represents.
#' @param famM : (datatype) Description of structure. What it represents.
#'
#' @return a list of two propensities
#'
get_propensity <- function( data ,
                            mX ,
                            mM ,
                            famX ,
                            famM ) {
  X<-glm(mX,data=data,family=famX)$fitted.values
  M<-glm(mM,data=data,family=famM)$fitted.values
  propensity <- list(X, M)
  return(propensity)
}

#  get_vars
#' This is an an internal function called \code{get_vars()}
#'
#' This returns a dataframe from vectors for exposure, mediator, c_xy, c_my, and
#' outcome as well as the individual vectors for exposure, mediator, and
#' outcome. It also checks to see if there are NULL values forany of these and
#' stops the program if any values are NULL.
#'
#' @keywords internal
#'
#' @param x : (datatype) Description of structure. What it represents.
#' @param m : (datatype) Description of structure. What it represents.
#' @param c_xy : (datatype) Description of structure. What it represents.
#' @param c_my : (datatype) Description of structure. What it represents.
#' @param y : (datatype) Description of structure. What it represents.
#' @param mY : (datatype) Description of structure. What it represents.
#' @param mM : (datatype) Description of structure. What it represents.
#' @param mX : (datatype) Description of structure. What it represents.
#'
#' @return a list composed of one data frame and three vectors
#'
get_vars <- function(x, m, c_xy, c_my, y, mX, mM, mY) {

  if(is.null(x)|is.null(m)|is.null(c_xy)|is.null(c_my)|is.null(y)){
    stop("Vectors are required for the exposure, mediator, c_xy, c_my, and outcome")
  }

  if(is.null(mX)|is.null(mM)|is.null(mY)){
    stop("Models are required for the outcome, exposure, and mediator")
  }

  df <- data.frame(x, m, c_xy, c_my, y)
  names(df) <- c("x", "m", "c_xy", "c_my", "y")

  vars <- list(df, x$x, m$m, y$y)

  return(vars)
}

#' Controlled Direct Effect (CDE) of the inverse probability-
#' weighted (IPW) marginal structural model.
#'
#' @param x : (datatype) Description of structure. What it represents.
#' @param m : (datatype) Description of structure. What it represents.
#' @param c_xy : (datatype) Description of structure. What it represents.
#' @param c_my : (datatype) Description of structure. What it represents.
#' @param y : (datatype) Description of structure. What it represents.
#' @param mX : (datatype) Description of structure. What it represents.
#' @param mM : (datatype) Description of structure. What it represents.
#' @param mY : (datatype) Description of structure. What it represents.
#' @param famX : (datatype) Description of structure. What it represents.
#' @param famM : (datatype) Description of structure. What it represents.
#' @param famY : (datatype) Description of structure. What it represents.
#'
#' @return : (vector) a length-2 list with the CDE with RR-interaction and
#'   RD-interaction.
#' @import zeallot
#' @references Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M.
#'   (2016). Mediation Analysis for Health Disparities Research. American
#'   Journal of Epidemiology, 184(4), 315–324.
#'   https://doi.org/10.1093/aje/kwv329
#' @author Ashley Naimi \email{ashley.naimi@pitt.edu}
#' @author Cantwell Carson \email{carsonc@gmail.com}
#' @examples 'examples/cde_example.R'
cde_ipw<-function(x=NULL,
                  m=NULL,
                  c_xy=NULL,
                  c_my=NULL,
                  y=NULL,
                  mX=NULL,
                  mM=NULL,
                  mY=NULL,
                  famX="gaussian",
                  famM="gaussian",
                  famY="gaussian"){

  c(data, exposure, mediator, outcome) %<-% get_vars(x=x, m=m, c_xy=c_xy, c_my=c_my, y=y,
                                                     mX=mX, mM=mM, mY=mY)

  # potential for superlearner...
  c(X_propensity, M_propensity) %<-% get_propensity( data=data,
                                                     mX=mX,
                                                     mM=mM,
                                                     famX=famX,
                                                     famM=famM )

  den <- (X_propensity * exposure + (1 - X_propensity) * (1 - exposure)) *
    (M_propensity * mediator + (1 - M_propensity) * (1 - mediator))

  num <- (mean(exposure) * exposure + (1 - mean(exposure)) * (1 - exposure)) *
    (mean(mediator) * mediator + (1 - mean(mediator)) * (1 - mediator))

  data$sw = num / den

  RR <- exp(coef(glm(mY, data=data, family=poisson(link="log"),
                     weights=sw))[as.character(mX[[2]])])

  RD <- coef(glm(mY, data=data, family=gaussian(link="identity"),
                 weights=sw))[as.character(mX[[2]])]

  ipw <- cbind(RD,RR)

  names(RR)<-"risk ratio"
  names(RD) <- "risk difference"
  row.names(ipw)<-"IPW"

  return(ipw)

}

#' Controlled Direct Effect (CDE) of the structural
#' transformation model.
#'
#' @param x : (datatype) Description of structure. What it represents.
#' @param m : (datatype) Description of structure. What it represents.
#' @param c_xy : (datatype) Description of structure. What it represents.
#' @param c_my : (datatype) Description of structure. What it represents.
#' @param y : (datatype) Description of structure. What it represents.
#' @param mX : (datatype) Description of structure. What it represents.
#' @param mM : (datatype) Description of structure. What it represents.
#' @param mY : (datatype) Description of structure. What it represents.
#' @param famX : (datatype) Description of structure. What it represents.
#' @param famM : (datatype) Description of structure. What it represents.
#' @param famY : (datatype) Description of structure. What it represents.
#'
#' @return : (vector) a length-2 list with the CDE with RR-interaction and
#'   RD-interaction.
#' @import zeallot
#' @references Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M.
#'   (2016). Mediation Analysis for Health Disparities Research. American
#'   Journal of Epidemiology, 184(4), 315–324.
#'   https://doi.org/10.1093/aje/kwv329
#' @author Ashley Naimi \email{ashley.naimi@pitt.edu}
#' @author Cantwell Carson \email{carsonc@gmail.com}
#' @examples 'examples/cde_example.R'
cde_transf<-function(x=NULL,
                     m=NULL,
                     c_xy=NULL,
                     c_my=NULL,
                     y=NULL,
                     mX=NULL,
                     mM=NULL,
                     mY=NULL,
                     famX="gaussian",
                     famM="gaussian",
                     famY="gaussian"){

  c(data, exposure, mediator, outcome) %<-% get_vars(x=x, m=m, c_xy=c_xy, c_my=c_my, y=y,
                                                     mX=mX, mM=mM, mY=mY)

  mod1 <- matrix(coef(lm(mY, data=data))[c("m","x:m")])
  y_tilde <- outcome - matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1
  RD <- coef(lm(y_tilde~exposure+c_xy,data=data))["exposure"]

  mod1 <- matrix(coef(glm(mY,data=data,family=poisson))[c("m","x:m")])
  y_tilde <- outcome*exp(- matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1)
  RR <- exp(coef(glm(y_tilde~exposure+c_xy,data=data,family=quasipoisson))["exposure"])

  transf<-cbind(RD,RR)

  names(RD)<-"Risk Difference"
  names(RR)<-"Risk Ratio"
  row.names(transf)<-"Struct Transf"

  return(transf)

}

#' Controlled Direct Effect (CDE) of the G-estimation model.

#' @param x : (datatype) Description of structure. What it represents.
#' @param m : (datatype) Description of structure. What it represents.
#' @param c_xy : (datatype) Description of structure. What it represents.
#' @param c_my : (datatype) Description of structure. What it represents.
#' @param y : (datatype) Description of structure. What it represents.
#' @param mX : (datatype) Description of structure. What it represents.
#' @param mM : (datatype) Description of structure. What it represents.
#' @param mY : (datatype) Description of structure. What it represents.
#' @param famX : (datatype) Description of structure. What it represents.
#' @param famM : (datatype) Description of structure. What it represents.
#' @param famY : (datatype) Description of structure. What it represents.
#'
#' @return : (vector) a length-2 list with the CDE with RR-interaction and
#'   RD-interaction.
#' @import zeallot
#' @references Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M.
#'   (2016). Mediation Analysis for Health Disparities Research. American
#'   Journal of Epidemiology, 184(4), 315–324.
#'   https://doi.org/10.1093/aje/kwv329
#' @author Ashley Naimi \email{ashley.naimi@pitt.edu}
#' @author Cantwell Carson \email{carsonc@gmail.com}
#' @examples 'examples/cde_example.R'
cde_gest<-function(x=NULL,
                   m=NULL,
                   c_xy=NULL,
                   c_my=NULL,
                   y=NULL,
                   mX=NULL,
                   mM=NULL,
                   mY=NULL,
                   famX="gaussian",
                   famM="gaussian",
                   famY="gaussian"){

  c(data, exposure, mediator, outcome) %<-% get_vars(x=x, m=m, c_xy=c_xy, c_my=c_my, y=y,
                                                     mX=mX, mM=mM, mY=mY)

  c(X_propensity, M_propensity) %<-% get_propensity(data=data, mX=mM, mM=mM, famX=famX, famM=famM )

  mod1<-matrix(coef(lm(y~x+I(m-M_propensity)+x:I(m-M_propensity)+c_my,data=data))[c("I(m - M_propensity)","x:I(m - M_propensity)")])
  y_tilde<-outcome - matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1
  RD<-coef(lm(y_tilde~I(x-X_propensity)+c_xy,data=data))["I(x - X_propensity)"]

  mod1<-matrix(coef(glm(y~x+I(m-M_propensity)+x:I(m-M_propensity)+c_xy+c_my,family=poisson,data=data))[c("I(m - M_propensity)","x:I(m - M_propensity)")])
  y_tilde<-outcome*exp( - matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1)
  RR<-exp(coef(glm(y_tilde~I(x-X_propensity)+c_xy,family=quasipoisson,data=data))["I(x - X_propensity)"])

  gest<-cbind(RD,RR)

  names(RD)<-"Risk Difference"
  names(RR)<-"Risk Ratio"
  row.names(gest)<-"G Estim"

  return(gest)
}

#' Controlled Direct Effect (CDE) of the Targeted Minimum Loss-based
#' Estimation (TMLE) model.
#'
#' @param x : (datatype) Description of structure. What it represents.
#' @param m : (datatype) Description of structure. What it represents.
#' @param c_xy : (datatype) Description of structure. What it represents.
#' @param c_my : (datatype) Description of structure. What it represents.
#' @param y : (datatype) Description of structure. What it represents.
#' @param mX : (datatype) Description of structure. What it represents.
#' @param mM : (datatype) Description of structure. What it represents.
#' @param mY : (datatype) Description of structure. What it represents.
#' @param famX : (datatype) Description of structure. What it represents.
#' @param famM : (datatype) Description of structure. What it represents.
#' @param famY : (datatype) Description of structure. What it represents.
#'
#' @return : (vector) a length-2 list with the CDE with RR-interaction and
#'   RD-interaction.
#' @import zeallot
#' @references Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M.
#'   (2016). Mediation Analysis for Health Disparities Research. American
#'   Journal of Epidemiology, 184(4), 315–324.
#'   https://doi.org/10.1093/aje/kwv329
#' @author Ashley Naimi \email{ashley.naimi@pitt.edu}
#' @author Cantwell Carson \email{carsonc@gmail.com}
#' @examples 'examples/cde_example.R'
cde_tmle<-function(x=NULL,
                   m=NULL,
                   c_xy=NULL,
                   c_my=NULL,
                   y=NULL,
                   mX=NULL,
                   mM=NULL,
                   mY=NULL,
                   famX="gaussian",
                   famM="gaussian",
                   famY="gaussian"){

  c(data, exposure, mediator, outcome) %<-% get_vars(x=x, m=m, c_xy=c_xy, c_my=c_my, y=y,
                                                     mX=mX, mM=mM, mY=mY)

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
  tml<-cbind(RD,RR)

  names(RD)<-"Risk Difference"
  names(RR)<-"Risk Ratio"
  row.names(tml)<-"TMLE"

  return(tml)

}

#' Controlled Direct Effect (CDE) model comparison
#'
#' @param x : (datatype) Description of structure. What it represents.
#' @param m : (datatype) Description of structure. What it represents.
#' @param c_xy : (datatype) Description of structure. What it represents.
#' @param c_my : (datatype) Description of structure. What it represents.
#' @param y : (datatype) Description of structure. What it represents.
#' @param mX : (datatype) Description of structure. What it represents.
#' @param mM : (datatype) Description of structure. What it represents.
#' @param mY : (datatype) Description of structure. What it represents.
#' @param famX : (datatype) Description of structure. What it represents.
#' @param famM : (datatype) Description of structure. What it represents.
#' @param famY : (datatype) Description of structure. What it represents.
#'
#' @return : (data frame) a (2, 4) data frame with the RR-interaction and
#'   RD-interaction for each of the models included in the analysis.
#' @references Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M.
#'   (2016). Mediation Analysis for Health Disparities Research. American
#'   Journal of Epidemiology, 184(4), 315–324.
#'   https://doi.org/10.1093/aje/kwv329
#' @author Ashley Naimi \email{ashley.naimi@pitt.edu}
#' @author Cantwell Carson \email{carsonc@gmail.com}
#' @examples 'examples/cde_example.R'
compare_cde<-function(x=NULL,
                      m=NULL,
                      c_xy=NULL,
                      c_my=NULL,
                      y=NULL,
                      mX=NULL,
                      mM=NULL,
                      mY=NULL,
                      famX="gaussian",
                      famM="gaussian",
                      famY="gaussian"){

  args = list(x, m, c_xy, c_my, y, mX, mM, mY, famX, famM, famY)

  estimators = list(cde_ipw, cde_transf, cde_gest, cde_tmle)

  cde_vals <- data.frame()

  for (estimator in estimators){
    cde_vals <- rbind(cde_vals, do.call(estimator, args))
  }

  return(cde_vals)

}
