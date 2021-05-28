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



#datalong format
datalong = list()
for (i in 1:N.sim){
  dat1 <- data1[[i]] %>%select(studlab,r1,n1,t1)
  colnames(dat1) <- c("study", "responders", "sampleSize", "treatment")
  
  dat2 <- data1[[i]] %>%select(studlab,r2,n2,t2)
  colnames(dat2) <- c("study", "responders", "sampleSize", "treatment")
  
  datalong[[i]] <- rbind(dat1, dat2)
  
}


for (i in 1:N.sim){ datalong[[i]]=list("data"=datalong[[i]],"logOR"=logOR[[i]]) }



IV=function(X)
{
  biasIV.RE = c()
  coverageIV.RE = c()
  TAU2.RE = c()
  d = c()
  CId = c()
  mtc.network <- mtc.network(data.ab = X$data, description = "Network")
  mtc.model1 <- mtc.model(mtc.network, type ="consistency", likelihood = "binom", link = "logit", 
                          hy.prior=mtc.hy.prior("prec", "dgamma", 0.001, 0.001),  n.chain = "2",linearModel = "random", om.scale = 5)
  results1 <- mtc.run(mtc.model1, n.adapt=7500, n.iter=22500)
  x =summary(results1)$summaries
  Bay.res = data.frame(mean = x$statistics[1:(NT-1),1], lowerCI = x$quantiles[1:(NT-1),1], upperCI =x$quantiles[1:(NT-1),5])
  biasIV.RE<-c(biasIV.RE, (Bay.res$mean-X$logOR))
  Bay.res$cover<-(Bay.res$lowerCI<X$logOR)&(Bay.res$upperCI>X$logOR)
  coverageIV.RE<-c(coverageIV.RE, Bay.res$cover)
  TAU2.RE <- c(TAU2.RE,x$statistics[NT]^2)
  for (i in 1:(NT-1)) {
    d[i] = c(Bay.res$mean[i])}
  for (i in 1:(NT-1)) {
    CId[i] = c(Bay.res$upperCI[i]-Bay.res$lowerCI[i])}
  return(list("biasRE"=biasIV.RE,"covRE"=coverageIV.RE, "Tau2RE"=TAU2.RE, "d" = d,"CId" = CId))
}


# Initiate cluster
no_cores <- detectCores()
cl <- makeCluster(no_cores)
clusterExport(cl, "IV")
clusterExport(cl,"datalong")
clusterExport(cl,"NT")
clusterExport(cl,"N.sim")
clusterExport(cl, "IV")
clusterEvalQ(cl, {library(gemtc)})
clusterEvalQ(cl, {library(rjags)})
l1=parLapply(cl,1:N.sim, function(x) IV(datalong[[x]]))


#Treatment Estimate
d = c()
for (i in 1:N.sim){  d=c(d, l1[[i]]$d)}
d_NMA_BayG = c()
for (i in 1:(NT-1)) {
  d_NMA_BayG[i] =mean(d[seq(i, length(d), NT-1)])}
D_NMA_BayG[r] = list(d_NMA_BayG)

#Confidence Interval
CId = c()
for (i in 1:N.sim){  CId=c(CId, l1[[i]]$CId)}
CId_NMA_BayG = c()
for (i in 1:(NT-1)) {
  CId_NMA_BayG[i] =mean(CId[seq(i, length(CId), NT-1)])}
CID_NMA_BayG[r] = list(CId_NMA_BayG)
MeanCI_NMA_BayG[r] = mean(CId)

#Bias and RMSE
biasIV.RE=c()
for (i in 1:N.sim){  biasIV.RE=c(biasIV.RE, l1[[i]]$biasRE)}
BIAS_NMA_BayG[r] = mean(biasIV.RE) #Bias
RMSE_NMA_BayG[r]= sqrt(mean((biasIV.RE)^2)) #RMSE

#Coverage probability
coverageIV.RE=c()
for (i in 1:N.sim){  coverageIV.RE=c(coverageIV.RE, l1[[i]]$covRE)}
CP_NMA_BayG[r] = mean(coverageIV.RE) #CP

#Tau squared
TAU2.RE=c()
for (i in 1:N.sim){  TAU2.RE=c(TAU2.RE, l1[[i]]$Tau2RE)}
TAU2_NMA_BayG[r] = mean(TAU2.RE) #TAU2

#End cluster
stopCluster(cl)
