library(trimediation)
library(mediation)

data("framing", package="mediation")
example_data <- framing
set.seed(2014)

mX <- x~c_xy
mM <- m~x+c_xy+c_my
mY <- y~x+m+x:m+c_xy+c_my

x    <- example_data["tone"]
m    <- example_data["eth"]
c_xy <- example_data["anti_info"]
c_my <- example_data["treat"]
y    <- example_data["cong_mesg"]

colnames(x)    <- "x"
colnames(m)    <- "m"
colnames(c_xy) <- "c_xy"
colnames(c_my) <- "c_my"
colnames(y)    <- "y"

## cde_ipw and cde_tmle both work with this model. cde_transf and cde_gest do
## not. I don't know enough to fix them. Also, I don't know how to handle the
## warnings in these instances.

cde_ipw(x=x,
            m=m,
            c_xy=c_xy,
            c_my=c_my,
            y=y,
            mX=mX,
            mM=mM,
            mY=mY,
            famX="binomial",
            famM="binomial",
            famY="binomial",
            boot=TRUE,
            sims = 100)
