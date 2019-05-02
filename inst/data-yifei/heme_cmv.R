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
