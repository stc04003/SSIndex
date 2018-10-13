library(GSM)

dim(simDat(100, "M1"))
dim(simDat(100, "M2"))
dim(simDat(100, "M3"))
dim(simDat(100, "M4"))

c(7, 24) / 25

do <- function(n, model) {
    dat <- simDat(n, model)
    unlist(gsm(dat))
}

do(200, "M2")

sim1 <- t(replicate(100, do(200, "M1")))
sim2 <- t(replicate(100, do(200, "M2")))
sim3 <- t(replicate(100, do(200, "M3")))
sim4 <- t(replicate(100, do(200, "M4")))

summary(sim1)
##      b01               b02               r01              r02        
## Min.   :-0.9996   Min.   :-1.0000   Min.   :0.1063   Min.   :0.9140  
## 1st Qu.:-0.4778   1st Qu.:-0.4956   1st Qu.:0.2424   1st Qu.:0.9501  
## Median : 0.3819   Median : 0.4404   Median :0.2760   Median :0.9611  
## Mean   : 0.1879   Mean   : 0.2056   Mean   :0.2738   Mean   :0.9601  
## 3rd Qu.: 0.7694   3rd Qu.: 0.8214   3rd Qu.:0.3121   3rd Qu.:0.9702  
## Max.   : 1.0000   Max.   : 1.0000   Max.   :0.4056   Max.   :0.9943  

summary(sim2)
##       b01               b02               r01              r02        
##  Min.   :-0.9992   Min.   :-0.9998   Min.   :0.3832   Min.   :0.6246  
##  1st Qu.:-0.6802   1st Qu.:-0.8523   1st Qu.:0.5701   1st Qu.:0.7566  
##  Median :-0.6073   Median :-0.7945   Median :0.6143   Median :0.7891  
##  Mean   :-0.5882   Mean   :-0.7549   Mean   :0.6139   Mean   :0.7852  
##  3rd Qu.:-0.5231   3rd Qu.:-0.7331   3rd Qu.:0.6538   3rd Qu.:0.8215  
##  Max.   : 0.4551   Max.   : 0.4454   Max.   :0.7809   Max.   :0.9237 

summary(sim3)
##       b01               b02               r01               r02         
##  Min.   :-0.9984   Min.   :-0.9991   Min.   :-0.7811   Min.   :-0.9328  
##  1st Qu.:-0.6610   1st Qu.:-0.8424   1st Qu.:-0.6500   1st Qu.:-0.8314  
##  Median :-0.6075   Median :-0.7943   Median :-0.6040   Median :-0.7970  
##  Mean   :-0.5849   Mean   :-0.7709   Mean   :-0.6011   Mean   :-0.7938  
##  3rd Qu.:-0.5389   3rd Qu.:-0.7504   3rd Qu.:-0.5556   3rd Qu.:-0.7599  
##  Max.   : 0.3865   Max.   : 0.4133   Max.   :-0.3605   Max.   :-0.6244  

summary(sim4)
##       b01               b02               r01                r02         
##  Min.   :-0.9996   Min.   :-1.0000   Min.   :-1.00000   Min.   :-0.2680  
##  1st Qu.:-0.6576   1st Qu.:-0.8425   1st Qu.:-0.87818   1st Qu.: 0.4783  
##  Median :-0.5974   Median :-0.8019   Median :-0.76504   Median : 0.6440  
##  Mean   :-0.5877   Mean   :-0.7745   Mean   :-0.71497   Mean   : 0.6197  
##  3rd Qu.:-0.5386   3rd Qu.:-0.7534   3rd Qu.:-0.60696   3rd Qu.: 0.7947  
##  Max.   : 0.2988   Max.   : 0.3511   Max.   : 0.05221   Max.   : 0.9999  

Douglas06(simDat(100, "M2")) ## confirms data generation


do(200, "M1")
do(200, "M2")
do(200, "M3")
do(200, "M4")
