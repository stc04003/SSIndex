library(GSM)
library(reReg)
library(frailtyHL)

data(cgd)
## ID 87 ends with status = 1
cgd <- rbind(cgd, cgd[cgd$id == 87,][2,])
cgd$status[nrow(cgd)] <- 0
cgd <- cgd[order(cgd$id),]
rownames(cgd) <- NULL

with(cgd, reSurv(tstop, id, status))

str(gsm(reSurv(tstop, id, status) ~ propylac, data = cgd))
