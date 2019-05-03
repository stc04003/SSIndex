## setwd("~/Dropbox/Marr/data")

all <- read.csv("all.ser.inf.long.time.dep.hsct.dat.csv",header=T) 
all$I.heme_state_not_remission<-rep(1,dim(all)[1])
all$I.heme_state_not_remission[all$heme_state==1]<-0 ### remission
all$I.heme.state.relapse<-rep(NA,dim(all)[1])
all$I.heme.state.relapse[all$heme_state!=2]<-0
all$I.heme.state.relapse[all$heme_state==2]<-1
## 0 D+R+
## 1 D-R+
## 2 D+R-
## 3 D-R-
## For (D+/R-, D-/R-) vs (D+/R+, D-/R+)
all$I.cmv2<-rep(NA,dim(all)[1])
all$I.cmv2[all$I.cmv.d.r>=2]<-0
all$I.cmv2[all$I.cmv.d.r<2]<-1
all$I.cmv2[all$I.cmv.d.r==4]<-NA
## For D-/R- vs (D+/R+. D-/R+, D+/R-)
all$I.cmv.dneg.rneg<-rep(NA,dim(all)[1])
all$I.cmv.dneg.rneg[all$I.cmv.d.r!=3]<-1
all$I.cmv.dneg.rneg[all$I.cmv.d.r==3]<-0
all$I.cmv.dneg.rneg[all$I.cmv.d.r==4]<-NA

all

all$study_id

################################################################################################
## Take 2
################################################################################################
all <- read.csv("hsct.long_2014_09_24.csv",header=T)
all$heme_state

length(unique(all$study_id)) ## 174
sum(!is.na(all$heme_state)) ## 151
sum(!is.na(all$sero_cmv_hsct)) ## 171
sum(!is.na(all$dsero_cmv_hsct)) ## 153

table(all$heme_state)
table(all$sero_cmv_hsct)
table(all$dsero_cmv_hsct)

all$study_id[!is.na(all$heme_state)]
all$study_id[!is.na(all$sero_cmv_hsct)]
all$study_id[!is.na(all$dsero_cmv_hsct)]

table(all$sero_cmv_hsct, all$dsero_cmv_hsct)


heme <- data.frame(id = all$study_id[!is.na(all$heme_state)],
                   heme = all$heme_state[!is.na(all$heme_state)])

sero <- data.frame(id = all$study_id[!is.na(all$sero_cmv_hsct)],
                   sero = all$sero_cmv_hsct[!is.na(all$sero_cmv_hsct)])

dsero <- data.frame(id = all$study_id[!is.na(all$dsero_cmv_hsct)],
                    dsero = all$dsero_cmv_hsct[!is.na(all$dsero_cmv_hsct)])

dim(sero)
head(sero)

## save(heme, file = "heme.RData")
## save(sero, file = "sero.RData")
## save(dsero, file = "dsero.RData")
