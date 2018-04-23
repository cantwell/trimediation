##########################
# This script is intended to provide a demonstration of the controlled direct
# effect (CDE) as applied to four different models as with random
# data. The four different models are:
#    1) Inverse Probability-Weighted (IPW) marginal structural model
#    2) Structural Transformation model
#    3) G-estimation model
#    4) Targeted Minimum Loss-based Estimation (TMLE) model
#
# When run, it should produce the following:
#                        RD        RR
# IPW           -0.0004150921 1.087889
# Struct Transf  0.0486749783 1.631896
# G Estim        0.0266387378 1.359153
# TMLE           0.0464612669 1.540893
#
#
# These are described in:
#   Naimi, A. I., Schnitzer, M. E., Moodie, E. E. M., & Bodnar, L. M. (2016).
#   Mediation Analysis for Health Disparities Research. American Journal of
#   Epidemiology, 184(4), 315–324. https://doi.org/10.1093/aje/kwv329
# For more information, please contact ashley.naimi@pitt.edu
##########################

library(trimediation)

data(trimed_example)
example_data <-trimed_example
set.seed(2014)

mX <- x~c_xy
mM <- m~x+c_xy+c_my
mY <- y~x+m+x:m+c_xy+c_my

compare_cde(x=trimed_example["x"],
            m=trimed_example["m"],
            c_xy=trimed_example["c_xy"],
            c_my=trimed_example["c_my"],
            y=trimed_example["y"],
            mX=mX,
            mM=mM,
            mY=mY,
            famX="binomial",
            famM="binomial",
            famY="binomial",
            boot=TRUE,
            sims = 50)
