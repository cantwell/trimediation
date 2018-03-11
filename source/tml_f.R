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
#' @seealso \code{\link{cde_f}}
#' @seealso \code{\link{ipw_f}}
#' @seealso \code{\link{stran_f}}
#' @seealso \code{\link{g_est_f}}
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
