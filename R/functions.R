## 1) revised last functions using cde_ipw as template -- DONE
## 2) push everything to GH --DONE
## 3) implement bootstrap for CI --DONE
## 4) test functons in applied datasets --IN PROCESS

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
get_propensity <- function( tcde_data ,
                            tmX ,
                            tmM ,
                            tfamX ,
                            tfamM ) {
  X<-glm(tmX,data=tcde_data,family=tfamX)$fitted.values
  M<-glm(tmM,data=tcde_data,family=tfamM)$fitted.values
  propensity <- list(X, M)
  return(propensity)
}

#  get_resample
#' This is an an internal function called \code{get_resample()}
#'
#' This returns the rows of a data frame resampled with repetition for
#' bootstrapping boot intervals.
#'
#' @keywords internal
#'
#' @param df : (dataframe) An arbitrary dataframe to be resampled.
#'
#' @return a list of two propensities
#'
get_resample <- function(df){
  df <- df[sample(nrow(df), nrow(df), replace = TRUE),]
  return(df)
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

#  get_boot
#' This is an an internal function called \code{get_boot()}
#'
#' This returns a boot interval from a data set using manual bootstrap resampling.
#'
#' @keywords internal
#'
#' @param cde_ci : a dataframe containing the data necessary
#' @param ci : (float) A number between zero and 1 denoting the boot
#'   interval desired
#'
#' @return a number corresponding to the normal boot half-interval for the
#'   resampled dataset
#'
get_boot <- function(cde_ci, ci) {


  cde_mean <- apply(cde_ci, 2, mean)
  cde_std  <- apply(cde_ci, 2, sd)
  cde_len  <- nrow(cde_ci)
  cde_err <- qnorm(ci)*cde_std/sqrt(cde_len)

  return(cde_err)
}

#  get_output
#' This is an an internal function called \code{get_output()}
#'
#' This formats a controlled direct effect estimate for risk difference and
#' risk ratio along with the boot intervals from resampling (if requested).
#'
#' @keywords internal
#'
#' @param boot : (boolean) whether or not the boot interval should be calculated
#' @param df : a dataframe containing the data necessary

#' @param ci : (float) A number between zero and 1 denoting the boot
#'   interval desired
#'
#' @return a cbind list of risk ratio and risk difference along with their boot intervals
#'
get_output <- function(cde_ci, ci) {

  RR <- cde_ci[1,1]
  RD <- cde_ci[1,2]
  if (nrow(cde_ci > 1)) {

    cde_err <- get_boot(cde_ci, ci)

    RR_l_error <- RR - cde_err[1]
    RR_r_error <- RR + cde_err[1]
    RD_l_error <- RD - cde_err[2]
    RD_r_error <- RD + cde_err[2]

    cde_out <- cbind(RR, RR_l_error, RR_r_error, RD, RD_l_error, RD_r_error)
    names(RR)<-"risk ratio"
    names(RD) <- "risk difference"
    names(RD_l_error)<-"RD lower ci"
    names(RD_r_error) <- "RD upper ci"
    names(RR_l_error)<-"RR lower ci"
    names(RR_r_error) <- "RR upper ci"
  }

  else{
    cde_out <- cbind(RR, RD)
    names(RR)<-"risk ratio"
    names(RD) <- "risk difference"
  }

  return(cde_out)
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
#' @param boot : (boolean) A boolean variable that indicates whether a
#'   non-parametric bootstrap will be used to calculate confidence intervals
#' @param sims : (positive integer) If boot is TRUE, sims is the number of
#'   simulations used to generate the confidence intervals.
#' @param ci   : (float) A float between 0 and 1 representing the desired
#'   confidence interval to be used for the non-parametric bootstrap.
#'
#' @return : (vector) a length-2 list with the CDE with RR-interaction and
#'   RD-interaction. If boot is TRUE, it will also return the right and left
#'   bounds of the desired confidence interval, ci.
#'
#' @import zeallot
#' @export
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
                  famY="gaussian",
                  boot=TRUE,
                  sims=100,
                  ci=0.95){

  c(cde_data_o, exposure_o, mediator_o, outcome_o) %<-% get_vars(x=x, m=m, c_xy=c_xy, c_my=c_my, y=y,
                                                     mX=mX, mM=mM, mY=mY)
  cde_ci <- data.frame()

  if (!boot) {
    sims <- 1
  }

  for (i in 1:sims){

    if (i == 1){
      cde_data <- cde_data_o
      exposure <- exposure_o
      mediator <- mediator_o
      outcome  <- outcome_o
    } else {
      cde_data <- get_resample(cde_data_o)
      exposure <- cde_data$x
      x        <- cde_data$x
      mediator <- cde_data$m
      m        <- cde_data$m
      outcome  <- cde_data$y
      y        <- cde_data$y
      c_xy <- cde_data$c_xy
      c_my <- cde_data$c_my
    }

    # potential for superlearner...
    c(X_propensity, M_propensity) %<-% get_propensity( cde_data,
                                                       mX,
                                                       mM,
                                                       famX,
                                                       famM )

    den <- (X_propensity * exposure + (1 - X_propensity) * (1 - exposure)) *
      (M_propensity * mediator + (1 - M_propensity) * (1 - mediator))

    num <- (mean(exposure) * exposure + (1 - mean(exposure)) * (1 - exposure)) *
      (mean(mediator) * mediator + (1 - mean(mediator)) * (1 - mediator))

    cde_data$sw = num / den

    RR <- exp(coef(glm(mY, data=cde_data, family=poisson(link="log"),
                       weights=sw))[as.character(mX[[2]])])

    RD <- coef(glm(mY, data=cde_data, family=gaussian(link="identity"),
                   weights=sw))[as.character(mX[[2]])]

    cde_ci <- rbind(cde_ci, c(RR, RD))
  }

  output <- get_output(cde_ci, ci)
  row.names(output) <- "IPW"

  return(output)

}
###TODO test get_output and format to other functions.
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
#' @param boot : (boolean) A boolean variable that indicates whether a
#'   non-parametric bootstrap will be used to calculate confidence intervals
#' @param sims : (positive integer) If boot is TRUE, sims is the number of
#'   simulations used to generate the confidence intervals.
#' @param ci   : (float) A float between 0 and 1 representing the desired
#'   confidence interval to be used for the non-parametric bootstrap.
#'
#' @return : (vector) a length-2 list with the CDE with RR-interaction and
#'   RD-interaction. If boot is TRUE, it will also return the right and left
#'   bounds of the desired confidence interval, ci.
#'
#' @import zeallot
#' @export
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
                     famY="gaussian",
                     boot=TRUE,
                     sims=100,
                     ci=0.95){

  c(cde_data_o, exposure_o, mediator_o, outcome_o) %<-% get_vars(x=x, m=m, c_xy=c_xy, c_my=c_my, y=y,
                                                                 mX=mX, mM=mM, mY=mY)

  cde_ci <- data.frame()

  if (!boot) {
    sims <- 1
  }

  for (i in 1:sims){

    if (i == 1){
      cde_data <- cde_data_o
      exposure <- exposure_o
      mediator <- mediator_o
      outcome  <- outcome_o
    } else {
      cde_data <- get_resample(cde_data_o)
      exposure <- cde_data$x
      x        <- cde_data$x
      mediator <- cde_data$m
      m        <- cde_data$m
      outcome  <- cde_data$y
      y        <- cde_data$y
      c_xy <- cde_data$c_xy
      c_my <- cde_data$c_my

    }


    mod1 <- matrix(coef(lm(mY, data=cde_data))[c("m","x:m")])
    y_term = matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1
    y_tilde <- outcome - y_term
    RD <- coef(lm(y_tilde~exposure+c_xy,data=cde_data))["exposure"]

    mod1 <- matrix(coef(glm(mY,data=cde_data,family=poisson))[c("m","x:m")])
    y_term = matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1
    y_tilde <- outcome*exp(-y_term)
    RR <- exp(coef(glm(y_tilde~exposure+c_xy,data=cde_data,family=quasipoisson))["exposure"])

    cde_ci <- rbind(cde_ci, c(RR, RD))
  }

  output <- get_output(cde_ci, ci)
  row.names(output) <- "TRANSF"
  return(output)

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
#' @param boot : (boolean) A boolean variable that indicates whether a
#'   non-parametric bootstrap will be used to calculate confidence intervals
#' @param sims : (positive integer) If boot is TRUE, sims is the number of
#'   simulations used to generate the confidence intervals.
#' @param ci   : (float) A float between 0 and 1 representing the desired
#'   confidence interval to be used for the non-parametric bootstrap.
#'
#' @return : (vector) a length-2 list with the CDE with RR-interaction and
#'   RD-interaction. If boot is TRUE, it will also return the right and left
#'   bounds of the desired confidence interval, ci.
#'
#' @import zeallot
#' @export
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
                   famY="gaussian",
                   boot=TRUE,
                   sims=100,
                   ci=0.95){

  c(cde_data_o, exposure_o, mediator_o, outcome_o) %<-% get_vars(x=x, m=m, c_xy=c_xy, c_my=c_my, y=y,
                                                                 mX=mX, mM=mM, mY=mY)

  cde_ci <- data.frame()

  if (!boot) {
    sims <- 1
  }

  for (i in 1:sims){

    if (i == 1){
      cde_data <- cde_data_o
      exposure <- exposure_o
      mediator <- mediator_o
      outcome  <- outcome_o
    } else {
      cde_data <- get_resample(cde_data_o)
      exposure <- cde_data$x
      x        <- cde_data$x
      mediator <- cde_data$m
      m        <- cde_data$m
      outcome  <- cde_data$y
      y        <- cde_data$y
      c_xy <- cde_data$c_xy
      c_my <- cde_data$c_my

    }

    c(X_propensity, M_propensity) %<-% get_propensity( cde_data,
                                                       mX,
                                                       mM,
                                                       famX,
                                                       famM )

    mod1<-matrix(coef(lm(y~x+I(m-M_propensity)+x:I(m-M_propensity)+c_my,data=cde_data))[c("I(m - M_propensity)","x:I(m - M_propensity)")])
    y_tilde<-outcome - matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1
    RD<-coef(lm(y_tilde~I(x-X_propensity)+c_xy,data=cde_data))["I(x - X_propensity)"]

    mod1<-matrix(coef(glm(y~x+I(m-M_propensity)+x:I(m-M_propensity)+c_xy+c_my,family=poisson,data=cde_data))[c("I(m - M_propensity)","x:I(m - M_propensity)")])
    y_tilde<-outcome*exp( - matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1)
    RR<-exp(coef(glm(y_tilde~I(x-X_propensity)+c_xy,family=quasipoisson,data=cde_data))["I(x - X_propensity)"])

    cde_ci <- rbind(cde_ci, c(RR, RD))
  }

  output <- get_output(cde_ci, ci)
  row.names(output) <- "GEST"
  return(output)
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
#' @param boot : (boolean) A boolean variable that indicates whether a
#'   non-parametric bootstrap will be used to calculate confidence intervals
#' @param sims : (positive integer) If boot is TRUE, sims is the number of
#'   simulations used to generate the confidence intervals.
#' @param ci   : (float) A float between 0 and 1 representing the desired
#'   confidence interval to be used for the non-parametric bootstrap.
#'
#' @return : (vector) a length-2 list with the CDE with RR-interaction and
#'   RD-interaction. If boot is TRUE, it will also return the right and left
#'   bounds of the desired confidence interval, ci.
#'
#' @import zeallot
#' @export
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
                   famY="gaussian",
                   boot=TRUE,
                   sims=100,
                   ci=0.95){

  c(cde_data_o, exposure_o, mediator_o, outcome_o) %<-% get_vars(x=x, m=m, c_xy=c_xy, c_my=c_my, y=y,
                                                                 mX=mX, mM=mM, mY=mY)

  cde_ci <- data.frame()

  if (!boot) {
    sims <- 1
  }

  for (i in 1:sims){

    if (i == 1){
      cde_data <- cde_data_o
      exposure <- exposure_o
      mediator <- mediator_o
      outcome  <- outcome_o
    } else {
      cde_data <- get_resample(cde_data_o)
      exposure <- cde_data$x
      x        <- cde_data$x
      mediator <- cde_data$m
      m        <- cde_data$m
      outcome  <- cde_data$y
      y        <- cde_data$y
      c_xy <- cde_data$c_xy
      c_my <- cde_data$c_my

    }

    mY<-glm(y~c_xy+c_my,data=cde_data,subset=x==1&m==0,family=binomial(link="logit"))
    mM<-glm(m~c_my,data=cde_data,subset=x==1,family=binomial(link="logit"))
    mX<-glm(x~c_xy,data=cde_data,family=binomial(link="logit"))

    newdata1<-cde_data;newdata1$x<-1;newdata1$m<-0;
    pY<-predict(mY,newdata=newdata1,type="response")
    pM<-predict(mM,newdata=newdata1,type="response")
    pX<-predict(mX,newdata=newdata1,type="response")
    cde_data$cc2 <- as.numeric(cde_data$x==1&cde_data$m==0)/((1-pM)*(pX))
    q2_star<-glm(y~-1+cc2,data=cde_data,offset=log(pY/(1-pY)),family=binomial(link="logit"))$fitted.values
    mQ1<-glm(q2_star ~ c_xy,data=cde_data,subset=x==1,family=quasibinomial(link="logit"))
    q1<-predict(mQ1,newdata=newdata1,type="response")
    cde_data$cc1 <- as.numeric(cde_data$x==1)/pX
    mu10<-glm(q2_star ~ -1 + cc1, data=cde_data,family=quasibinomial(link="logit"),offset=log(q1/(1-q1)))$fitted.values

    mY<-glm(y~c_xy+c_my,data=cde_data,subset=x==0&m==0,family=binomial(link="logit"))
    mM<-glm(m~c_my,data=cde_data,subset=x==0,family=binomial(link="logit"))

    newdata1<-cde_data;newdata1$x<-0;newdata1$m<-0;
    pY<-predict(mY,newdata=newdata1,type="response")
    pM<-predict(mM,newdata=newdata1,type="response")
    pX<-predict(mX,newdata=newdata1,type="response")
    cde_data$cc2 <- as.numeric(cde_data$x==0&cde_data$m==0)/((1-pM)*(pX))
    q2_star<-glm(y~-1+cc2,data=cde_data,offset=log(pY/(1-pY)),family=binomial(link="logit"))$fitted.values
    mQ1<-glm(q2_star ~ c_xy,data=cde_data,subset=x==0,family=quasibinomial(link="logit"))
    q1<-predict(mQ1,newdata=newdata1,type="response")
    cde_data$cc1 <- as.numeric(cde_data$x==0)/pX
    mu00<-glm(q2_star ~ -1 + cc1, data=cde_data,family=quasibinomial(link="logit"),offset=log(q1/(1-q1)))$fitted.values

    RD<-mean(mu10-mu00)
    RR<-mean(mu10)/mean(mu00)

    cde_ci <- rbind(cde_ci, c(RR, RD))
  }

  output <- get_output(cde_ci, ci)

  row.names(output) <- "TMLE"
  return(output)
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
#' @param boot : (boolean) A boolean variable that indicates whether a
#'   non-parametric bootstrap will be used to calculate confidence intervals
#' @param sims : (positive integer) If boot is TRUE, sims is the number of
#'   simulations used to generate the confidence intervals.
#' @param ci   : (float) A float between 0 and 1 representing the desired
#'   confidence interval to be used for the non-parametric bootstrap.
#' @export
#'
#' @return : (data frame) a data frame with all four mediation methods for
#'   compaisionm, with RR-interaction and RD-interaction. If boot is TRUE, it
#'   will also return the right and left bounds of the desired confidence
#'   interval, ci.
#'
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
                      famY="gaussian",
                      boot=TRUE,
                      sims=100,
                      ci=0.95){

  args = list(x, m, c_xy, c_my, y, mX, mM, mY, famX, famM, famY, boot, sims, ci)

  estimators = list(cde_ipw, cde_transf, cde_gest, cde_tmle)

  cde_vals <- data.frame()

  for (estimator in estimators){
    cde_vals <- rbind(cde_vals, do.call(estimator, args))
  }

  return(cde_vals)

}
