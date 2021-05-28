rm(list=ls())
set.seed(42) # the answer to the ultimate question of life, the universe and everything
N.sim=10
library(MASS)


data1
attach(data1[[10]])
data1[[10]]$name = c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10)
library(metafor)
rma(odds.t1, p.t1, data=data1[[10]],  mods=~factor(name),method="REML")

detach(data1[[10]])
x

name="Scenario-19.RData"
data1=list()
logOR=list()
logOR1=list()
OR=list()

NT=5 #### number of treatments in the network
NS=2 #### number of studies per comparison
tau=0.1 ####  heterogeneity SD
Npmin=100 #### minimum number of patients per arm
Npmax=200 #### maximum number of patients per arm

### define treatment indices
t1=c()
t2=c()
for (i in 1:(NT-1)){
  for (k in (i+1):NT){
    for(j in 1:NS){
      t1=c(t1,i)
      t2=c(t2,k)      }}}
N.stud=length(t1)

### define patients per treatment arm
for (i in 1:N.sim)
{   
  logOR[[i]]=seq(from =1/(NT-1), to = 1, by = 1/(NT-1))  ### true logOR across studies  ### true logOR across studies 
  OR[[i]]=c(1,exp(logOR[[i]]))
  data1[[i]]=data.frame(t1,t2)
  data1[[i]]$studlab=c(1:(N.stud))
  data1[[i]]$n1=data1[[i]]$n2=round(runif(N.stud,Npmin,Npmax))}


#### define probabilities per treatment, per study arm
for (i in 1:N.sim)
{  
  data1[[i]]$p.ref=runif(N.stud,0.3,0.5) #### study-specific probability of an event in treatment 1
  data1[[i]]$odds.ref=data1[[i]]$p.ref/(1-data1[[i]]$p.ref)
}

#### define probabilities per treatment, per study arm
Sigma=matrix(c(tau^2,tau^2/2,tau^2/2,tau^2),nrow=2)

for (i in 1:N.sim)
{   
  logOR1[[i]]=c(0,logOR[[i]])
  for(j in 1:N.stud){
    
    
    data1[[i]]$truelogOR.t1[j]=logOR1[[i]][data1[[i]]$t1[j]]
    data1[[i]]$truelogOR.t2[j]=logOR1[[i]][data1[[i]]$t2[j]]
    test1=mvrnorm(1,c(data1[[i]]$truelogOR.t1[j],data1[[i]]$truelogOR.t2[j]),Sigma)
    test2=rnorm(1,data1[[i]]$truelogOR.t2[j],tau)
    data1[[i]]$st.sp.logOR.t1[j]=test1[1]*(data1[[i]]$t1[j]!=1)
    data1[[i]]$st.sp.logOR.t2[j]=test1[2]*(data1[[i]]$t1[j]!=1)+test2*(data1[[i]]$t1[j]==1)
    
    data1[[i]]$odds.t1[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t1[j])
    data1[[i]]$odds.t2[j]=data1[[i]]$odds.ref[j]*exp(data1[[i]]$st.sp.logOR.t2[j])
    data1[[i]]$p.t1[j]=data1[[i]]$odds.t1[j]/(1+data1[[i]]$odds.t1[j])
    data1[[i]]$p.t2[j]=data1[[i]]$odds.t2[j]/(1+data1[[i]]$odds.t2[j])
  }}


#### generate the data
for (i in 1:N.sim)
{  
  for(j in 1:N.stud){
    data1[[i]]$r1[j]=rbinom(1,data1[[i]]$n1[j],data1[[i]]$p.t1[j]) 
    data1[[i]]$r2[j]=rbinom(1,data1[[i]]$n2[j],data1[[i]]$p.t2[j]) 
  }}

for (i in 1:N.sim){ data1[[i]]=data1[[i]][,-c(6:11)]}
##################

x = analyseNetmeta(data1[[2]])
x

analyseNetmeta <- function(data){
  
  # calculate the log-OR and its standard error
  data <- data %>%
    rowwise() %>%
    mutate(TE = log((r2/(n2 - r2))/
                      (r1/(n1 - r1))), 
           seTE = sqrt((1/r2) + (1/r1) + 
                         (1/(n2 - r2)) + (1/(n1 - r1)))) 
  
  # perform the network meta-analysis 
  netmeta.model <- netmeta(data = data, TE = TE, seTE = seTE, 
                           treat1 = t1, treat2 = t2, 
                           studlab = studlab, sm = "OR")
  
  return(netmeta.model)
}

netgraph(x)





X1=list()
for (i in 1:N.sim){ X1[[i]]=list("data"=data1[[i]],"logOR"=logOR[[i]]) }
X1





IV=function(X)
{
  biasIV.FE<-c()
  coverageIV.FE<-c()
  biasIV.RE<-c()
  coverageIV.RE<-c()
  
  IV1<-netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "Inverse", incr=0.5 ,cc.pooled=T, allstudies=T ,sm="OR",data=X$data)
  IV.FE.res<-data.frame(mean=IV1$TE.fixed[2:NT,1], lowerCI=IV1$lower.fixed[2:NT,1],upperCI=IV1$upper.fixed[2:NT,1])
  biasIV.FE<-c(biasIV.FE, (IV.FE.res$mean-X$logOR))
  IV.FE.res$cover<-(IV.FE.res$lowerCI<X$logOR)&(IV.FE.res$upperCI>X$logOR)
  coverageIV.FE<-c(coverageIV.FE, IV.FE.res$cover)
  
  IV.RE.res<-data.frame(mean=IV1$TE.random[2:NT,1], lowerCI=IV1$lower.random[2:NT,1],upperCI=IV1$upper.random[2:NT,1])
  biasIV.RE<-c(biasIV.RE, (IV.RE.res$mean-X$logOR))
  IV.RE.res$cover<-(IV.RE.res$lowerCI<X$logOR)&(IV.RE.res$upperCI>X$logOR)
  coverageIV.RE<-c(coverageIV.RE, IV.RE.res$cover)
  
  return(list("biasFE"=biasIV.FE,"biasRE"=biasIV.RE,"covFE"=coverageIV.FE, "covRE"=coverageIV.RE))
  
}
library(parallel)
install.packages('beepr')
library(beepr)


# Initiate cluster
no_cores <- detectCores()
cl <- makeCluster(no_cores)
clusterExport(cl, "IV")
clusterExport(cl,"X1")
clusterExport(cl,"NT")
clusterExport(cl,"N.sim")
clusterExport(cl, "IV")
clusterEvalQ(cl, {library(netmeta)})
l1=parLapply(cl,1:N.sim, function(x) IV(X1[[x]]))


biasIV.FE=c()
for (i in 1:N.sim){  biasIV.FE=c(biasIV.FE, l1[[i]]$biasFE)}
mean(biasIV.FE)

coverageIV.FE=c()
for (i in 1:N.sim){  coverageIV.FE=c(coverageIV.FE, l1[[i]]$covFE)}
mean(coverageIV.FE)

biasIV.RE=c()
for (i in 1:N.sim){  biasIV.RE=c(biasIV.RE, l1[[i]]$biasRE)}
mean(biasIV.RE)

coverageIV.RE=c()
for (i in 1:N.sim){  coverageIV.RE=c(coverageIV.RE, l1[[i]]$covRE)}
mean(coverageIV.RE)

stopCluster(cl)
beep(sound=3)
