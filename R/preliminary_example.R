
d<-data.frame(x=c(rep(0,16),rep(1,16)),
              m=c(rep(0,8),rep(1,8),rep(0,8),rep(1,8)),
              c_xy=rep(c(rep(0,4),rep(1,4)),4),
              c_my=rep(c(0,0,1,1),8),
              y=rep(c(0,1),16),
              freq=c(707,21,223,31,442,23,122,39,309,
                     19,206,45,178,23,136,48,147,8,
                     192,23,216,15,248,76,176,15,371,
                     85,217,14,441,184))

D <- d[rep(row.names(d), d$freq),1:5];row.names(D)<-NULL

head(D);table(D$y)

mX<-as.formula(x~c_xy)
mM<-as.formula(m~x+c_xy+c_my)
mY<-as.formula(y~x+m+x:m+c_xy+c_my)

outcome<-D$y
mediator<-D$m
exposure<-D$x

famX=binomial
famM=binomial

dat=D

cde<-function(outcome,mY,famY,
              mediator,mM,famM,
              exposure,mX,famX,
              data,interactionRR=T,interactionRD=T){

  # IPW
  X_propensity<-glm(mX,data=dat,family=famX)$fitted.values
  M_propensity<-glm(mM,data=dat,family=famM)$fitted.values

  den<-(X_propensity*exposure + (1-X_propensity)*(1-exposure))*(M_propensity*mediator + (1-M_propensity)*(1-mediator))
  num<-(mean(exposure)*exposure + (1-mean(exposure))*(1-exposure))*(mean(mediator)*mediator + (1-mean(mediator))*(1-mediator))

  sw=num/den

  if(interactionRD==T){
    RD<-coef(glm(outcome~exposure+mediator+exposure:mediator,data=dat,family=gaussian(link="identity"),weights=sw))["exposure"]
    names(RD)<-"risk difference"
  } else{
    RD<-coef(glm(outcome~exposure+mediator,data=dat,family=gaussian(link="identity"),weights=sw))["exposure"]
    names(RD)<-"risk difference"
  }

  if(interactionRR==T){
    RR<-exp(coef(glm(outcome~exposure+mediator+exposure:mediator,data=dat,family=poisson(link="log"),weights=sw))["exposure"])
    names(RR)<-"risk ratio"
  } else{
    RR<-exp(coef(glm(outcome~exposure+mediator,data=dat,family=poisson(link="log"),weights=sw))["exposure"])
    names(RR)<-"risk ratio"
  }

  ipw<-cbind(RD,RR);row.names(ipw)<-"IPW"

  # structural transformation
  if(interactionRD==T){
    mod1<-matrix(coef(lm(y~x+m+x:m+c_xy+c_my,data=dat))[c("m","x:m")])
    y_tilde<-outcome - matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1
    RD<-coef(lm(y_tilde~exposure+c_xy,data=dat))["exposure"]
  } else{
    mod1<-matrix(coef(lm(y~x+m+c_xy+c_my,data=dat))[c("m")])
    y_tilde<-outcome - matrix(mediator,ncol=1)%*%mod1
    RD<-coef(lm(y_tilde~exposure+c_xy,data=dat))["exposure"]
  }

  if(interactionRR==T){
    mod1<-matrix(coef(glm(y~x+m+x:m+c_xy+c_my,data=dat,family=poisson))[c("m","x:m")])
    y_tilde<-outcome*exp(- matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1)
    RR<-exp(coef(glm(y_tilde~exposure+c_xy,data=dat,family=quasipoisson))["exposure"])
  } else{
    mod1<-matrix(coef(glm(y~x+m+c_xy+c_my,data=dat,family=poisson))[c("m")])
    y_tilde<-outcome*exp( - matrix(mediator,ncol=1)%*%mod1)
    RR<-exp(coef(glm(y_tilde~exposure+c_xy,data=dat,family=poisson))["exposure"])
  }

  names(RD)<-"Risk Difference"
  names(RR)<-"Risk Ratio"
  stran<-cbind(RD,RR);row.names(stran)<-"Struct Transf"

  # G - estimation
  if(interactionRD==T){
    mod1<-matrix(coef(lm(y~x+I(m-M_propensity)+x:I(m-M_propensity)+c_my,data=dat))[c("I(m - M_propensity)","x:I(m - M_propensity)")])
    y_tilde<-outcome - matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1
    RD<-coef(lm(y_tilde~I(x-X_propensity)+c_xy,data=dat))["I(x - X_propensity)"]
  } else{
    mod1<-matrix(coef(lm(y~x+I(m-M_propensity)+c_my,data=dat))[c("I(m - M_propensity)")])
    y_tilde<-outcome - matrix(mediator,ncol=1)%*%mod1
    RD<-coef(lm(y_tilde~I(x-X_propensity)+c_xy,data=dat))["I(x - X_propensity)"]
  }

  if(interactionRR==T){
    mod1<-matrix(coef(glm(y~x+I(m-M_propensity)+x:I(m-M_propensity)+c_xy+c_my,family=poisson,data=dat))[c("I(m - M_propensity)","x:I(m - M_propensity)")])
    y_tilde<-outcome*exp( - matrix(c(mediator,mediator*exposure),ncol=2)%*%mod1)
    RR<-exp(coef(glm(y_tilde~I(x-X_propensity)+c_xy,family=quasipoisson,data=dat))["I(x - X_propensity)"])
  } else{
    mod1<-matrix(coef(glm(y~x+I(m-M_propensity)+c_xy+c_my,family=poisson,data=dat))[c("I(m - M_propensity)")])
    y_tilde<-outcome*exp( - matrix(mediator,ncol=1)%*%mod1)
    RR<-exp(coef(lm(y_tilde~I(x-X_propensity)+c_xy,family=quasipoisson,data=dat))["I(x - X_propensity)"])
  }

  names(RD)<-"Risk Difference"
  names(RR)<-"Risk Ratio"
  g_est<-cbind(RD,RR);row.names(g_est)<-"G Estim"


  ## TMLE
  mY<-glm(y~c_xy+c_my,data=dat,subset=x==1&m==0,family=binomial(link="logit"))
  mM<-glm(m~c_my,data=dat,subset=x==1,family=binomial(link="logit"))
  mX<-glm(x~c_xy,data=dat,family=binomial(link="logit"))
  newDat1<-dat;newDat1$x<-1;newDat1$m<-0;
  pY<-predict(mY,newdata=newDat1,type="response")
  pM<-predict(mM,newdata=newDat1,type="response")
  pX<-predict(mX,newdata=newDat1,type="response")
  dat$cc2 <- as.numeric(dat$x==1&dat$m==0)/((1-pM)*(pX))
  q2_star<-glm(y~-1+cc2,data=dat,offset=log(pY/(1-pY)),family=binomial(link="logit"))$fitted.values
  mQ1<-glm(q2_star ~ c_xy,data=dat,subset=x==1,family=quasibinomial(link="logit"))
  q1<-predict(mQ1,newdata=newDat1,type="response")
  dat$cc1 <- as.numeric(dat$x==1)/pX
  mu10<-glm(q2_star ~ -1 + cc1, data=dat,family=quasibinomial(link="logit"),offset=log(q1/(1-q1)))$fitted.values

  mY<-glm(y~c_xy+c_my,data=dat,subset=x==0&m==0,family=binomial(link="logit"))
  mM<-glm(m~c_my,data=dat,subset=x==0,family=binomial(link="logit"))
  mX<-glm(x~c_xy,data=dat,family=binomial(link="logit"))
  newDat1<-dat;newDat1$x<-0;newDat1$m<-0;
  pY<-predict(mY,newdata=newDat1,type="response")
  pM<-predict(mM,newdata=newDat1,type="response")
  pX<-predict(mX,newdata=newDat1,type="response")
  dat$cc2 <- as.numeric(dat$x==0&dat$m==0)/((1-pM)*(pX))
  q2_star<-glm(y~-1+cc2,data=dat,offset=log(pY/(1-pY)),family=binomial(link="logit"))$fitted.values
  mQ1<-glm(q2_star ~ c_xy,data=dat,subset=x==0,family=quasibinomial(link="logit"))
  q1<-predict(mQ1,newdata=newDat1,type="response")
  dat$cc1 <- as.numeric(dat$x==0)/pX
  mu00<-glm(q2_star ~ -1 + cc1, data=dat,family=quasibinomial(link="logit"),offset=log(q1/(1-q1)))$fitted.values

  RD<-mean(mu10-mu00)
  RR<-mean(mu10)/mean(mu00)
  names(RD)<-"Risk Difference"
  names(RR)<-"Risk Ratio"
  tml<-cbind(RD,RR);row.names(tml)<-"TMLE"


  return(rbind(ipw,stran,g_est,tml))

}

cde(outcome=D$y,mY=as.formula(y~x+m+x:m+c_xy+c_my),binomial,
    mediator=D$m,mM=as.formula(m~x+c_xy+c_my),binomial,
    exposure=D$x,mX=as.formula(x~c_xy),binomial,
    dat=D,interactionRD=T,interactionRR=T)
