data1=list()
logOR=list()
logOR1=list()
OR=list()


Npmin=75 #### minimum number of patients per arm
Npmax=125 #### maximum number of patients per arm

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
  data1[[i]]$p.ref=runif(N.stud,0.2,0.3) #### study-specific probability of an event in treatment 1
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

#calc logOR and var
for (i in 1:N.sim){
  d <- data1[[i]] %>%
    mutate(LOGOR = log((r2/(n2 - r2))/
                         (r1/(n1 - r1))), 
           varTE = ((1/r2) + (1/r1) + 
                      (1/(n2 - r2)) + (1/(n1 - r1))))
  data1[[i]] = d
}


for (i in 1:N.sim) {
  data1[[i]]$name = paste("T",t1,t2, sep = "")
}




#Function
IV=function(X)
{
  biasIV.RE = c()
  coverageIV.RE = c()
  TAU2.RE = c()
  d = c()
  CId = c()
  LOGNT = function(x) sequence(x:1)/x
  logOR_rma = LOGNT(NT-1)
  RMA = rma(LOGOR, varTE , data=X,  mods=~-1+factor(name),method="REML")
  X.TE.res = coef(summary(RMA))[-c(2,3,4)]
  TAU2.RE <- c(TAU2.RE,RMA$tau2)
  biasIV.RE<-c(biasIV.RE, (X.TE.res$estimate-logOR_rma))
  X.TE.res$cover<-(X.TE.res$ci.lb<logOR_rma)&(X.TE.res$ci.ub>logOR_rma)
  coverageIV.RE<-c(coverageIV.RE, X.TE.res$cover)
  for (i in 1:(NT-1)) {
    d[i] = c(X.TE.res$estimate[i])}
  for (i in 1:(NT-1)) {
    CId[i] = c(X.TE.res$ci.ub[i]-X.TE.res$ci.lb[i])}
  return(list("biasRE"=biasIV.RE,"covRE"=coverageIV.RE, "Tau2RE"=TAU2.RE, "d" = d,"CId" = CId))
}



# Initiate cluster
no_cores <- detectCores()
cl <- makeCluster(no_cores)
clusterExport(cl, "IV")
clusterExport(cl,"data1")
clusterExport(cl,"NT")
clusterExport(cl,"N.sim")
clusterExport(cl, "IV")
clusterEvalQ(cl, {library(metafor)})
l1=parLapply(cl,1:N.sim, function(x) IV(data1[[x]]))


#Treatment Estimate
d = c()
for (i in 1:N.sim){  d=c(d, l1[[i]]$d)}
d_RMA = c()
for (i in 1:(NT-1)) {
  d_RMA[i] =mean(d[seq(i, length(d), NT-1)])}
D_RMA[r] = list(d_RMA)

#Confidence Interval
CId = c()
for (i in 1:N.sim){  CId=c(CId, l1[[i]]$CId)}
CId_RMA = c()
for (i in 1:(NT-1)) {
  CId_RMA[i] =mean(CId[seq(i, length(CId), NT-1)])}
CID_RMA[r] = list(CId_RMA)
MeanCI_RMA[r] = mean(CId)

#Bias and RMSE
biasIV.RE=c()
for (i in 1:N.sim){  biasIV.RE=c(biasIV.RE, l1[[i]]$biasRE)}
BIAS_RMA[r] = mean(biasIV.RE) #Bias
RMSE_RMA[r]= sqrt(mean((biasIV.RE)^2)) #RMSE

#Coverage probability
coverageIV.RE=c()
for (i in 1:N.sim){  coverageIV.RE=c(coverageIV.RE, l1[[i]]$covRE)}
CP_RMA[r] = mean(coverageIV.RE) #CP

#Tau squared
TAU2.RE=c()
for (i in 1:N.sim){  TAU2.RE=c(TAU2.RE, l1[[i]]$Tau2RE)}
TAU2_RMA[r] = mean(TAU2.RE) #TAU2

#End cluster
stopCluster(cl)
