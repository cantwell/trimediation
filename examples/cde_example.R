##########################
# This script is intended to provide a demonstration of the counter-factual
# disparity measure (CDM) as applied to four different models as with random
# data. The four different models are:
#    1) Inverse Probability-Weighted (IPW) marginal structural model
#    2) Structural Transformation model
#    3) G-estimation model
#    4) Targeted Minimum Loss-based Estimation (TMLE) model
# These are described in:
#   Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M. (2016).
#   Mediation Analysis for Health Disparities Research. American Journal of
#   Epidemiology, 184(4), 315–324. https://doi.org/10.1093/aje/kwv329
# For more information, please contact ashley.naimi@pitt.edu
##########################

library(trimediation)

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

famY=binomial
famX=binomial
famM=binomial

cde_f(outcome=D$y,
      mY=as.formula(y~x+m+x:m+c_xy+c_my),
      famY=binomial,
      mediator=D$m,
      mM=as.formula(m~x+c_xy+c_my),
      famM=binomial,
      exposure=D$x,
      mX=as.formula(x~c_xy),
      famX=binomial,
      data=D,
      interactionRD=T,
      interactionRR=T)
