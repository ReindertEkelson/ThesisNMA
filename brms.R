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

#make new label
for (i in 1:N.sim) {
  data1[[i]]$name = paste("T",t1,t2, sep = "")
}


IV=function(X)
{  
  biasIV.RE = c()
  coverageIV.RE = c()
  TAU2.RE = c()
  d = c()
  CId = c()
  LOGNT = function(x) sequence(x:1)/x
  logOR_rma = LOGNT(NT-1)
  logOR_rma_int = c(logOR_rma[1], logOR_rma[2:length(logOR_rma)] - 1/(NT-1))
  priors <- c(prior(normal(0,75^2), class = Intercept),
              prior(cauchy(0,0.5), class = sd))
  m.brm <- brm(LOGOR|se(sqrt(varTE)) ~  1 +  name + (1|studlab),
               data = X,
               prior = priors,
               iter = 30000, chains = 2, thin = 1, warmup = 7500, control = list(adapt_delta = 0.999))
  x = summary(m.brm)
  brm.TE.res = x$fixed[,c(1,3,4)]
  colnames(brm.TE.res) = cbind("Estimate", "lowerCI", "upperCI")
  brm.TE.res = as.data.frame(brm.TE.res)
  TAU2.RE = c(TAU2.RE, (x$random$studlab[1])^2)
  biasIV.RE<-c(biasIV.RE, (brm.TE.res$Estimate-logOR_rma_int))
  brm.TE.res$cover<-(brm.TE.res$lowerCI<logOR_rma_int)&(brm.TE.res$upperCI>logOR_rma_int)
  coverageIV.RE<-c(coverageIV.RE, brm.TE.res$cover)
  for (i in 1:(NT)) {
    d[i] = c(brm.TE.res$Estimate[i])}
  for (i in 1:(NT)) {
    CId[i] = c(brm.TE.res$upperCI[i]-brm.TE.res$lowerCI[i])}
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
clusterEvalQ(cl, {library(brms)})
l1=parLapply(cl,1:N.sim, function(x) IV(data1[[x]]))


#Treatment Estimate
d = c()
for (i in 1:N.sim){  d=c(d, l1[[i]]$d)}
d_BRMS = c()
for (i in 1:(NT)) {
  d_BRMS[i] =mean(d[seq(i, length(d), NT)])}
D_BRMS[r] = list(d_BRMS)

#Confidence Interval
CId = c()
for (i in 1:N.sim){  CId=c(CId, l1[[i]]$CId)}
CId_BRMS = c()
for (i in 1:(NT)) {
  CId_BRMS[i] =mean(CId[seq(i, length(CId), NT)])}
CID_BRMS[r] = list(CId_BRMS)
MeanCI_BRMS[r] = mean(CId)

#Bias and RMSE
biasIV.RE=c()
for (i in 1:N.sim){  biasIV.RE=c(biasIV.RE, l1[[i]]$biasRE)}
BIAS_BRMS[r] = mean(biasIV.RE) #Bias
RMSE_BRMS[r]= sqrt(mean((biasIV.RE)^2)) #RMSE

#Coverage probability
coverageIV.RE=c()
for (i in 1:N.sim){  coverageIV.RE=c(coverageIV.RE, l1[[i]]$covRE)}
CP_BRMS[r] = mean(coverageIV.RE) #CP

#Tau squared
TAU2.RE=c()
for (i in 1:N.sim){  TAU2.RE=c(TAU2.RE, l1[[i]]$Tau2RE)}
TAU2_BRMS[r] = mean(TAU2.RE) #TAU2

#End cluster
stopCluster(cl)
