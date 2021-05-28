#base
rm(list=ls()) 
set.seed(42) 
N.sim=2000
N.Scen = 27

#libraries
library(mixmeta)
library(progress)
library(brms)
library(MASS)
library(tidyverse)
library(netmeta)
library(parallel)
library(metafor)
library(beepr)
library(tictoc)
library(gemtc)
library(pcnetmeta)
library(mvmeta)
library(rjags)
library(plyr)
library(Rglpk)


#Master
DATA = expand.grid(data.frame( c(3,5,8), c(3,4,5), c(0.1,sqrt(0.1),1)))
colnames(DATA) = c("NS", "NT", "tau")
DATA
##########################################################################
tic()
#metafor
BIAS_RMA = RMSE_RMA = CP_RMA = D_RMA = CID_RMA = TAU2_RMA = MeanCI_RMA = c()

pb <- progress_bar$new(
  format = "  Simulating [:bar] :percent in :elapsed",
  total = nrow(DATA), clear = FALSE, width= 60)

for (r in 1:nrow(DATA)){
  tryCatch({
  NT = DATA[r,]$NT
  NS = DATA[r,]$NS
  tau = DATA[r,]$tau
  source("metafor.R")
  pb$tick()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  }
beep(sound = 8)

DF_Results = cbind(BIAS_RMA, RMSE_RMA, CP_RMA, MeanCI_RMA, TAU2_RMA)
rownames(DF_Results) = paste("Scenario", 1:dim(DF_Results)[1])

D_RMA[D_RMA == "NULL"] <- NA
DF_Estimates = plyr::ldply(D_RMA, rbind)
rownames(DF_Estimates) = paste("Scenario", 1:dim(DF_Estimates)[1])
?length
CID_RMA[CID_RMA == "NULL"] <- NA
DF_CI = plyr::ldply(CID_RMA, rbind)
rownames(DF_CI) = paste("Scenario", 1:dim(DF_CI)[1])

write.csv(DF_Results,"~/Thesis/DF_Results_metafor.csv")
write.csv(DF_Estimates,"~/Thesis/DF_Estimates_metafor.csv")
write.csv(DF_CI,"~/Thesis/DF_CI_metafor.csv")




#netmeta
BIAS_NMA = RMSE_NMA = CP_NMA = D_NMA = CID_NMA = TAU2_NMA = MeanCI_NMA = c()

pb <- progress_bar$new(
  format = "  Simulating [:bar] :percent in :elapsed",
  total = nrow(DATA), clear = FALSE, width= 60)

for (r in 1:nrow(DATA)){
  tryCatch({
    NT = DATA[r,]$NT
    NS = DATA[r,]$NS
    tau = DATA[r,]$tau
    source("netmeta.R")
    pb$tick()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
beep(sound = 8)

DF_Results = cbind(BIAS_NMA, RMSE_NMA, CP_NMA, MeanCI_NMA, TAU2_NMA)
rownames(DF_Results) = paste("Scenario", 1:N.Scen)

D_NMA[D_NMA == "NULL"] <- NA
DF_Estimates = plyr::ldply(D_NMA, rbind)
rownames(DF_Estimates) = paste("Scenario", 1:N.Scen)

CID_NMA[CID_NMA == "NULL"] <- NA
DF_CI = plyr::ldply(CID_NMA, rbind)
rownames(DF_CI) = paste("Scenario", 1:N.Scen)

write.csv(DF_Results,"~/Thesis/DF_Results_netmeta.csv")
write.csv(DF_Estimates,"~/Thesis/DF_Estimates_netmeta.csv")
write.csv(DF_CI,"~/Thesis/DF_CI_netmeta.csv")




#gemtcU
BIAS_NMA_BayU = RMSE_NMA_BayU = CP_NMA_BayU = D_NMA_BayU = CID_NMA_BayU = TAU2_NMA_BayU = MeanCI_NMA_BayU =c()

pb <- progress_bar$new(
  format = "  Simulating [:bar] :percent in :elapsed",
  total = nrow(DATA), clear = FALSE, width= 60)

for (r in 1:nrow(DATA)){
  tryCatch({
    NT = DATA[r,]$NT
    NS = DATA[r,]$NS
    tau = DATA[r,]$tau
    source("gemtcU.R")
    pb$tick()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
beep(sound = 8)

DF_Results = cbind(BIAS_NMA_BayU, RMSE_NMA_BayU, CP_NMA_BayU, MeanCI_NMA_BayU, TAU2_NMA_BayU)
rownames(DF_Results) = paste("Scenario", 1:N.Scen)

D_NMA_BayU[D_NMA_BayU == "NULL"] <- NA
DF_Estimates = plyr::ldply(D_NMA_BayU, rbind)
rownames(DF_Estimates) = paste("Scenario", 1:N.Scen)

CID_NMA_BayU[CID_NMA_BayU == "NULL"] <- NA
DF_CI = plyr::ldply(CID_NMA_BayU, rbind)
rownames(DF_CI) = paste("Scenario", 1:N.Scen)

write.csv(DF_Results,"~/Thesis/DF_Results_gemtcU.csv")
write.csv(DF_Estimates,"~/Thesis/DF_Estimates_gemtcU.csv")
write.csv(DF_CI,"~/Thesis/DF_CI_gemtcU.csv")




#gemtcG
BIAS_NMA_BayG = RMSE_NMA_BayG = CP_NMA_BayG = D_NMA_BayG = CID_NMA_BayG = TAU2_NMA_BayG = MeanCI_NMA_BayG =c()

pb <- progress_bar$new(
  format = "  Simulating [:bar] :percent in :elapsed",
  total = nrow(DATA), clear = FALSE, width= 60)

for (r in 1:nrow(DATA)){
  tryCatch({
    NT = DATA[r,]$NT
    NS = DATA[r,]$NS
    tau = DATA[r,]$tau
    source("gemtcG.R")
    pb$tick()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
beep(sound = 8)

DF_Results = cbind(BIAS_NMA_BayG, RMSE_NMA_BayG, CP_NMA_BayG, MeanCI_NMA_BayG, TAU2_NMA_BayG)
rownames(DF_Results) = paste("Scenario", 1:N.Scen)

D_NMA_BayG[D_NMA_BayG == "NULL"] <- NA
DF_Estimates = plyr::ldply(D_NMA_BayG, rbind)
rownames(DF_Estimates) = paste("Scenario", 1:N.Scen)

CID_NMA_BayG[CID_NMA_BayG == "NULL"] <- NA
DF_CI = plyr::ldply(CID_NMA_BayG, rbind)
rownames(DF_CI) = paste("Scenario", 1:N.Scen)

write.csv(DF_Results,"~/Thesis/DF_Results_gemtcG.csv")
write.csv(DF_Estimates,"~/Thesis/DF_Estimates_gemtcG.csv")
write.csv(DF_CI,"~/Thesis/DF_CI_gemtcG.csv")





#gemtcHN
BIAS_NMA_BayHN = RMSE_NMA_BayHN = CP_NMA_BayHN = D_NMA_BayHN = CID_NMA_BayHN = TAU2_NMA_BayHN = MeanCI_NMA_BayHN =c()

pb <- progress_bar$new(
  format = "  Simulating [:bar] :percent in :elapsed",
  total = nrow(DATA), clear = FALSE, width= 60)

for (r in 1:nrow(DATA)){
  tryCatch({
    NT = DATA[r,]$NT
    NS = DATA[r,]$NS
    tau = DATA[r,]$tau
    source("gemtcHN.R")
    pb$tick()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
beep(sound = 8)

DF_Results = cbind(BIAS_NMA_BayHN, RMSE_NMA_BayHN, CP_NMA_BayHN, MeanCI_NMA_BayHN, TAU2_NMA_BayHN)
rownames(DF_Results) = paste("Scenario", 1:N.Scen)

D_NMA_BayHN[D_NMA_BayHN == "NULL"] <- NA
DF_Estimates = plyr::ldply(D_NMA_BayHN, rbind)
rownames(DF_Estimates) = paste("Scenario", 1:N.Scen)

CID_NMA_BayHN[CID_NMA_BayHN == "NULL"] <- NA
DF_CI = plyr::ldply(CID_NMA_BayHN, rbind)
rownames(DF_CI) = paste("Scenario", 1:N.Scen)

write.csv(DF_Results,"~/Thesis/DF_Results_gemtcHN.csv")
write.csv(DF_Estimates,"~/Thesis/DF_Estimates_gemtcHN.csv")
write.csv(DF_CI,"~/Thesis/DF_CI_gemtcHN.csv")


#pcnetmeta
BIAS_PC = RMSE_PC = CP_PC = D_PC = CID_PC = MeanCI_PC = c()

pb <- progress_bar$new(
  format = "  Simulating [:bar] :percent in :elapsed",
  total = nrow(DATA), clear = FALSE, width= 60)

for (r in 1:nrow(DATA)){
  tryCatch({
    NT = DATA[r,]$NT
    NS = DATA[r,]$NS
    tau = DATA[r,]$tau
    source("pcnetmeta.R")
    pb$tick()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
beep(sound = 8)

DF_Results = cbind(BIAS_PC, RMSE_PC, CP_PC, MeanCI_PC)
rownames(DF_Results) = paste("Scenario", 1:N.Scen)

D_PC[D_PC == "NULL"] <- NA
DF_Estimates = plyr::ldply(D_PC, rbind)
rownames(DF_Estimates) = paste("Scenario", 1:N.Scen)

CID_PC[CID_PC == "NULL"] <- NA
DF_CI = plyr::ldply(CID_PC, rbind)
rownames(DF_CI) = paste("Scenario", 1:N.Scen)

write.csv(DF_Results,"~/Thesis/DF_Results_pcnetmeta.csv")
write.csv(DF_Estimates,"~/Thesis/DF_Estimates_pcnetmeta.csv")
write.csv(DF_CI,"~/Thesis/DF_CI_pcnetmeta.csv")





#mixmeta
BIAS_MVM = RMSE_MVM = CP_MVM = D_MVM = CID_MVM = MeanCI_MVM = TAU2_MVM = c()

pb <- progress_bar$new(
  format = "  Simulating [:bar] :percent in :elapsed",
  total = nrow(DATA), clear = FALSE, width= 60)

for (r in 1:nrow(DATA)){
  tryCatch({
    NT = DATA[r,]$NT
    NS = DATA[r,]$NS
    tau = DATA[r,]$tau
    source("mixmeta.R")
    pb$tick()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
beep(sound = 8)

DF_Results = cbind(BIAS_MVM, RMSE_MVM, CP_MVM, MeanCI_MVM, TAU2_MVM)
rownames(DF_Results) = paste("Scenario", 1:dim(DF_Results)[1])

D_MVM[D_MVM == "NULL"] <- NA
DF_Estimates = plyr::ldply(D_MVM, rbind)
rownames(DF_Estimates) = paste("Scenario", 1:dim(DF_Estimates)[1])

CID_MVM[CID_MVM == "NULL"] <- NA
DF_CI = plyr::ldply(CID_MVM, rbind)
rownames(DF_CI) = paste("Scenario", 1:dim(DF_CI)[1])

write.csv(DF_Results,"~/Thesis/DF_Results_mixmeta.csv")
write.csv(DF_Estimates,"~/Thesis/DF_Estimates_mixmeta.csv")
write.csv(DF_CI,"~/Thesis/DF_CI_mixmeta.csv")



#brms
BIAS_BRMS = RMSE_BRMS = CP_BRMS = D_BRMS = CID_BRMS = TAU2_BRMS = MeanCI_BRMS = c()

pb <- progress_bar$new(
  format = "  Simulating [:bar] :percent in :elapsed",
  total = nrow(DATA), clear = FALSE, width= 60)

for (r in 1:nrow(DATA)){
  tryCatch({
    NT = DATA[r,]$NT
    NS = DATA[r,]$NS
    tau = DATA[r,]$tau
    source("brms.R")
    pb$tick()
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
beep(sound = 8)

DF_Results = cbind(BIAS_BRMS, RMSE_BRMS, CP_BRMS, MeanCI_BRMS, TAU2_BRMS)
rownames(DF_Results) = paste("Scenario", 1:N.Scen)

D_BRMS[D_BRMS == "NULL"] <- NA
DF_Estimates = plyr::ldply(D_BRMS, rbind)
rownames(DF_Estimates) = paste("Scenario", 1:N.Scen)

CID_BRMS[CID_BRMS == "NULL"] <- NA
DF_CI = plyr::ldply(CID_BRMS, rbind)
rownames(DF_CI) = paste("Scenario", 1:N.Scen)

write.csv(DF_Results,"~/Thesis/DF_Results_brms.csv")
write.csv(DF_Estimates,"~/Thesis/DF_Estimates_brms.csv")
write.csv(DF_CI,"~/Thesis/DF_CI_brms.csv")



beep(sound = 7)
toc()