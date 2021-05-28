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

X1=list()
for (i in 1:N.sim){ X1[[i]]=list("data"=data1[[i]],"logOR"=logOR[[i]]) }

IV=function(X)
{
  biasIV.RE = c()
  coverageIV.RE = c()
  TAU2.RE = c()
  d = c()
  CId = c()
  IV1<-netmetabin(event1=r1,event2=r2,n1=n1,n2=n2, studlab = studlab,treat1=t1, treat2=t2, method = "Inverse", incr=0.5 ,cc.pooled=T, allstudies=T ,sm="OR",data=X$data)
  IV.RE.res<-data.frame(mean=IV1$TE.random[2:NT,1], lowerCI=IV1$lower.random[2:NT,1],upperCI=IV1$upper.random[2:NT,1])
  biasIV.RE<-c(biasIV.RE, (IV.RE.res$mean-X$logOR))
  IV.RE.res$cover<-(IV.RE.res$lowerCI<X$logOR)&(IV.RE.res$upperCI>X$logOR)
  coverageIV.RE<-c(coverageIV.RE, IV.RE.res$cover)
  TAU2.RE <- c(TAU2.RE,IV1$tau^2)
  for (i in 1:(NT-1)) {
    d[i] = c(IV.RE.res$mean[i])}
  for (i in 1:(NT-1)) {
    CId[i] = c(IV.RE.res$upperCI[i]-IV.RE.res$lowerCI[i])}
  return(list("biasRE"=biasIV.RE,"covRE"=coverageIV.RE, "Tau2RE"=TAU2.RE, "d" = d,"CId" = CId))
}



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


#Treatment Estimate
d = c()
for (i in 1:N.sim){  d=c(d, l1[[i]]$d)}
d_NMA = c()
for (i in 1:(NT-1)) {
  d_NMA[i] =mean(d[seq(i, length(d), NT-1)])}
D_NMA[r] = list(d_NMA)

#Confidence Interval
CId = c()
for (i in 1:N.sim){  CId=c(CId, l1[[i]]$CId)}
CId_NMA = c()
for (i in 1:(NT-1)) {
  CId_NMA[i] =mean(CId[seq(i, length(CId), NT-1)])}
CID_NMA[r] = list(CId_NMA)
MeanCI_NMA[r] = mean(CId)

#Bias and RMSE
biasIV.RE=c()
for (i in 1:N.sim){  biasIV.RE=c(biasIV.RE, l1[[i]]$biasRE)}
BIAS_NMA[r] = mean(biasIV.RE) #Bias
RMSE_NMA[r]= sqrt(mean((biasIV.RE)^2)) #RMSE

#Coverage probability
coverageIV.RE=c()
for (i in 1:N.sim){  coverageIV.RE=c(coverageIV.RE, l1[[i]]$covRE)}
CP_NMA[r] = mean(coverageIV.RE) #CP

#Tau squared
TAU2.RE=c()
for (i in 1:N.sim){  TAU2.RE=c(TAU2.RE, l1[[i]]$Tau2RE)}
TAU2_NMA[r] = mean(TAU2.RE) #TAU2

#End cluster
stopCluster(cl)

