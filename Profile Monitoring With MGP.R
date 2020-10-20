#' Profile Monitoring With MGP
#' Simulation
list.of.packages <- c("Matrix", "nloptr","microbenchmark")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(Matrix)
library(microbenchmark)
library(nloptr)

rm(list=ls())
########################################### Specify Input #################################################
n=2 #Number of outputs (measurement)
dp=30 #Number of Design Points -(1+x+x^2)
nsd=0.5 #Measurement Noise
fun=function(x,i) {
  if     (i==1){5+2*cos(x)}
  else if(i==2){-5+2*sin(1*x)}
  # else if(i==3){-10+exp((x%%2))*cos(2*x)}
  # else if(i==4){1+1*x^3}
}
########################################### Generate Data #################################################
trains=seq(0,2*pi,length.out=dp)
num=31
trainy1=lapply(1:n, function(i){lapply(1:num,function(j){fun(trains,i)+rnorm(dp,0,nsd)})}) #Initial Profiles
trainy2=matrix(, nrow=n,ncol=dp) #Mean matrix
for(i in 1:n){
  trainy2[i,]=apply(do.call("rbind",trainy1[[i]]),2,mean)
}
trainy=split(trainy2, row(trainy2)) #MGP Input#Staking the means# turn into a list
trainyn1=cbind(do.call('rbind',trainy1[[1]]),do.call('rbind',trainy1[[2]]))
trainyn=split(trainyn1, row(trainyn1))
y=unlist(trainy) #apply(trainyn1,2,mean)
leny=length(y)
######################################### Plot ###########################################
dev.off()
for(i in 1:num){
  plot(trains,trainy1[[1]][[i]],ylim=c(-10,10),xlim=c(0,2*pi),type="b",pch=2,lwd=2,xlab=NA,ylab=NA,cex.axis=1.2)
  par(new=T)
  }

##########################################################################################
################################## Kriging ###############################################
##########################################################################################
x=trains
num=31
##########################################  Covariance Quantitative ########################################
cyii=function(a,b,L)
{
  d=outer(a,b,`-`);I=outer(a,b,`==`)
  L[1]^2*exp(-0.5*d^2/L[2]^2)+ I*L[3]^2
}
##########################################  Full Covariance #################################################
C=function(x,H)
{
  cyii(x,x,H[c(1:3)])
}
##########################################  Mean Function ##################################################
mf=Vectorize(function(x) {0})
mu=mf(x)
######################################### Log Likelihood ####################################################
ypredk=list()
yvark=list()
######################################### Separate Kriging ##################################################
for(j in 1:n){
  ypredkin=list()
  yvarkin=list()
  for(i in 1:num){
    nn=length(trainy1[[j]][[i]])
    logL=function(H,fn)
    {
      B=C(x,H)
      deter=det(B)
      if(deter!=0) {a=0.5*(log(deter)+t(trainy1[[j]][[i]]-mu)%*%solve(B,trainy1[[j]][[i]]-mu)+log(2*pi)*nn)
      } else {
        ch=chol(B)
        logdeter=2*(sum(log(diag(ch))))
        a=0.5*(logdeter+t(trainy1[[j]][[i]]-mu)%*%solve(B,trainy1[[j]][[i]]-mu)+log(2*pi)*nn)
      }
      return(as.numeric(a))
    }
    logL_grad=function(H,fn)
    {
      return(nl.grad(H,fn))
    }
    ########################################### Optimize ########################################################
    x0=c(1,1,0.5)
    library(nloptr)
    opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 2000,print_level=3)
    t=proc.time()
    one=nloptr(x0=x0,eval_f=logL, eval_grad_f=logL_grad, opts= opts, fn=logL)
    proc.time()-t
    H0=one$solution
    
    ########################################## Prediction ######################################################
    xstar=x
    cyiiN=function(a,b,L)
    {
      d=outer(a,b,`-`);I=outer(a,b,`==`)
      L[1]^2*exp(-0.5*d^2/L[2]^2)
    }
    
    ## Prediction Mean
    Pk=cyiiN(xstar,x,H0[c(1:2)]) ## Close to C(xstar,x,H0)
    ypredkin[[i]]=mf(xstar)+Pk%*%solve(C(x,H0),trainy1[[j]][[i]]-mu)
    
    ## Prediction Variance
    Sk=cyii(xstar,xstar,H0[c(1:3)])
    
    yvarkin[[i]]=Sk  -  Pk%*%solve(C(x,H0),t(Pk))
  }
  ypredk[[j]]=apply(simplify2array(ypredkin), 1:2, mean)
  yvark[[j]]=apply(simplify2array(yvarkin), 1:2, mean)
}

########################################### Drawings ######################################################
dev.off()
par(mfrow=c(2,2))
for(i in 1:n){
  aa=seq(0,2*pi,length.out=50)
  plot(aa,fun(aa,i),xlim=c(0,2*pi),ylim=c(-10,10),type="l",pch=16,lwd=1,xlab=NA,ylab=NA,cex.axis=1.2,main=i)
  par(new=T)
  plot(xstar,ypredk[[i]],xlim=c(0,2*pi),ylim=c(-10,10),type="b",xlab=NA,ylab=NA,cex.axis=1.2,col='red')
  par(new=T)
  plot(xstar,ypredk[[i]]+3*sqrt(diag(yvark[[i]])),xlim=c(0,2*pi),ylim=c(-10,10),type="b",xlab=NA,ylab=NA,cex.axis=1.2,col='green')
  par(new=T)
  plot(xstar,ypredk[[i]]-3*sqrt(diag(yvark[[i]])),xlim=c(0,2*pi),ylim=c(-10,10),type="b",xlab=NA,ylab=NA,cex.axis=1.2,col='blue')
  points(trains,trainy[[i]],xlim=c(0,2*pi),ylim=c(-10,10),pch=2,lwd=2,xlab=NA,ylab=NA,cex.axis=1.2)
}
###########End of drawing #############

#Imporatnt Stuff: ypredk, yvark # which now are two lists containing different outputs
rm(H0,mu,nn,one,opts,x0,xstar,C,cyii,cyiiN,logL,logL_grad,mf,Pk,Sk,x,H00,H01)
###############################################################################################################
########################################### Multivaraite Gaussian Process #####################################
###############################################################################################################

##########################################Index Function##################################
index=function(n,m)
{
  p1=c() #row
  p2=c() #col
  for(i in 1:n){
    for(j in 1:m){
      p1=c(p1,(i-1)*m+seq(1:j))
      p2=c(p2,rep(j+(i-1)*m,j))
    }
  }
  p3=c()#row
  for(k in 1:(n-1)){
    for(i in (k+1):n){
      for(j in 1:m){
        p3=c(p3,(k-1)*m+seq(1:m))
      }
    }
  }
  
  p4=c() #col
  for(k in 1:(n-1)){
    for(i in (k+1):n){
      for(j in 1:m){
        p4=c(p4,rep(j+(i-1)*m,m))
      }
    }
  }
  return(list(pfi=c(p1,p3),pfj=c(p2,p4)))
}

pf=index(n,dp)
pfi=pf$pfi;pfj=pf$pfj 
# sparseMatrix(i=pfi,j=pfj) 
####################################### Covaraince Function #############################
cyii=function(a,b,L)
{
  d=outer(a,b,`-`);I=outer(a,b,`==`)
  d=d[upper.tri(d,diag=T)];I=I[upper.tri(I,diag=T)]
  L[1]^2*exp(-0.25*d^2/L[2]^2) + L[3]^2*exp(-0.25*d^2/L[4]^2) + I*L[5]^2
}
cyij=function(a,b,L)
{
  d=outer(a,b,`-`)
  L[1]*L[3]*sqrt(2*abs(L[2]*L[4])/(L[2]^2+L[4]^2))*exp(-0.5*d^2/(L[2]^2+L[4]^2))
}
###################################### Covariance Matrix #################################
C=function(strain,H){
  zii=list();zij=list()
  zii=lapply(1:n, function(i){cyii(strain,strain,H[c(2*i-1,2*i,2*n+2*i-1,2*n+2*i,4*n+1)])})
  zij=lapply(1:(n-1),function(i){lapply((i+1):n,function(j){cyij(strain,strain,H[c(2*n+2*i-1,2*n+2*i,2*n+2*j-1,2*n+2*j)])})})
  b1=unlist(zii);b2=unlist(zij)
  return(sparseMatrix(i=pfi,j=pfj,x=c(b1,b2),symmetric=T))
}
##################################### likelihood ###################################################(prod(diag(B)))^2***B=chol(C(s1,s2,H))***

ypredin=list()
yvarin=list()
for(i in 1:num){
  logL=function(H,fn)
  {
    B=C(trains,H)
    deter=det(B)
    if(deter>0) {a=0.5*(log(deter)+t(trainyn[[i]])%*%solve(B,trainyn[[i]])+log(2*pi)*leny)
    } else {
      ch=chol(B)
      logdeter=2*(sum(log(diag(ch))))
      a=0.5*(logdeter+t(trainyn[[i]])%*%solve(B,trainyn[[i]])+log(2*pi)*leny)
    }
    return(as.numeric(a))
  }
  logL_grad=function(H,fn)
  {
    return(nl.grad(H,fn))
  }
  
  ################################## Optimize #############################################
  x0=c(rep(1,4*n),rep(0.5,1))
  opts <- list( "algorithm" = "NLOPT_LD_MMA","maxeval" = 2500,print_level=3) #print_level=3
  #t=proc.time()
  one=tryCatch(nloptr(x0=x0,eval_f= logL,eval_grad_f = logL_grad,opts= opts,fn= logL ), error = function(e) e)
  if(any(class(one) == "error")==T){cat("nloptr Problem")}
  H1=one$solution
  five2=tryCatch(uobyqa(H1,fn=logL,control=list(maxfun=50000)),error = function(e) e) #iprint=3 in control=list..
  if(any(class(five2) == "error")==T){cat("uobyqa Problem")}
  H0=five2$par
  # result1=c(result1,five2$fval)
  
  ##################################Result#################################################
  # H0=apply(do.call("rbind",H01),2,mean,na.rm=TRUE)
  xstar=trains
  
  cyy=function(i,j,a,b,L,m){
    if(i==j){
      d=outer(a,b,`-`);I=outer(a,b,`==`)
      L[2*i-1]^2*exp(-0.25*d^2/L[2*i]^2) + L[2*n+2*i-1]^2*exp(-0.25*d^2/L[2*n+2*i]^2) + I*L[4*n+1]^2
    }
    else{
      d=outer(a,b,`-`);
      L[2*n+2*i-1]*L[2*n+2*j-1]*sqrt(2*abs(L[2*n+2*i]*L[2*n+2*j])/(L[2*n+2*i]^2+L[2*n+2*j]^2))*exp(-0.5*d^2/(L[2*n+2*i]^2+L[2*n+2*j]^2))
    }
  }
  
  sigma<-function(x,y,L,measures){
    d1=length(x);d2=length(y)
    sigma=matrix(,nrow=(measures*d1),ncol=(measures*d2))
    for(i in 1:measures){
      for(j in 1:measures){
        sigma[(((i-1)*d1+1):(i*d1)),(((j-1)*d2+1):(j*d2))]=cyy(i,j,x,y,L,m=measures)
      }
    }
    return(sigma)
  }
  #####New Sets of function for variance####
  cyyn=function(i,j,a,b,L,m){
    d=outer(a,b,`-`)
    if(i==j){
      L[2*i-1]^2*exp(-0.25*d^2/L[2*i]^2) + L[2*n+2*i-1]^2*exp(-0.25*d^2/L[2*n+2*i]^2)
    }
    else{
      L[2*n+2*i-1]*L[2*n+2*j-1]*sqrt(2*abs(L[2*n+2*i]*L[2*n+2*j])/(L[2*n+2*i]^2+L[2*n+2*j]^2))*exp(-0.5*d^2/(L[2*n+2*i]^2+L[2*n+2*j]^2))
    }
  }
  
  sigman<-function(x,y,L,measures){
    d1=length(x);d2=length(y)
    sigma=matrix(,nrow=(measures*d1),ncol=(measures*d2))
    for(i in 1:measures){
      for(j in 1:measures){
        sigma[(((i-1)*d1+1):(i*d1)),(((j-1)*d2+1):(j*d2))]=cyyn(i,j,x,y,L,m=measures)
      }
    }
    return(sigma)
  }
  
  etan=function(xs,x,L,measures){
    t(sigman(xs,x,L,measures))
  }
  ###########################################Predicting Mean and Variance ##################
  covM=C(trains,H0)
  pk=etan(xs=xstar,x=trains,L=H0,measures=n)
  sk=sigma(xstar,xstar,L=H0,measures=n)
  ypredin[[i]]=as.matrix(t(pk)%*%solve(covM,trainyn[[i]]))
  yvarin[[i]]=as.matrix(sk-t(pk)%*%solve(covM,pk))
  
}

ypred=apply(simplify2array(ypredin), 1:2, mean)
yvar=apply(simplify2array(yvarin), 1:2, mean)

lens=length(xstar)
ylist=list()
for(i in 1:n){
  ylist[[i]]=ypred[((i-1)*lens+1):(i*lens)]
}
varlist=list()
for(i in 1:n){
  varlist[[i]]=yvar[(((i-1)*lens+1):(i*lens)),(((i-1)*lens+1):(i*lens))]
}

for(i in 1:n){
  plot(aa,fun(aa,i),xlim=c(0,2*pi),ylim=c(-10,10),type="l",pch=16,lwd=1,xlab=NA,ylab=NA,cex.axis=1.2,main=i)
  par(new=T)
  plot(xstar,ylist[[i]],xlim=c(0,2*pi),ylim=c(-10,10),type="b",xlab=NA,ylab=NA,cex.axis=1.2,col='red')
  par(new=T)
  plot(xstar,ylist[[i]]+3*sqrt(diag(varlist[[i]])),xlim=c(0,2*pi),ylim=c(-10,10),type="b",xlab=NA,ylab=NA,cex.axis=1.2,col='green')
  par(new=T)
  plot(xstar,ylist[[i]]-3*sqrt(diag(varlist[[i]])),xlim=c(0,2*pi),ylim=c(-10,10),type="b",xlab=NA,ylab=NA,cex.axis=1.2,col='blue')
  points(trains,trainy[[i]],xlim=c(0,2*pi),ylim=c(-10,10),pch=2,lwd=2,xlab=NA,ylab=NA,cex.axis=1.2)
}

rm(etan,sigman,cyyn,sigma,cyy,logL_grad,logL,C,cyij,cyii)


#############################
#####Normal T-Sq Charting####
#############################
x=list()
for(i in 1:n){
  x[[i]]=do.call("rbind",trainy1[[i]])
}

sy=list()
for(i in 1:n){
  sy[[i]]=cov(x[[i]],x[[i]])
}
############################################################################################
############################### Profile Monitoring Simulation ##############################
############################################################################################
################################################
### Estimating Control Limits in both Charts ###
################################################
simulation=100000
t1=list()
t2=list()
t3=list()
m=1
while(m<=simulation){
  ##Generating new in-control profile ##
  yn=lapply(1:n, function(i){fun(trains,i)+rnorm(dp,0,nsd)})
  ##Forming T^2 Statistic##  Kriging  ###
  t1[[m]]=unlist(lapply(1:n, function(i){mu1=0;mu1=yn[[i]]-ypredk[[i]];t(mu1)%*%solve(yvark[[i]],mu1)}))
  
  ##Forming T^2 Statistic## MGP ###
  mu2=unlist(yn)-ypred
  t2[m]=as.numeric(t(mu2)%*%solve(yvar,mu2))
  
  ##Forming T^2 Statistic##Typical T^2###
  t3[[m]]=unlist(lapply(1:n, function(i){mu3=0;mu3=yn[[i]]-trainy[[i]];t(mu3)%*%solve(sy[[i]],mu3)}))
  print(m)
  m=m+1
}

q1=apply(do.call("rbind",t1),2,quantile,0.9973^(1/n))
q2=quantile(unlist(t2),0.9973)
q3=apply(do.call("rbind",t3),2,quantile,0.9973^(1/n))

##################
###Out of Control ARL Simulation###
# ARL1 for Kriging, ARL2 for MGP, ARL3 for T^2
gpm=list()
mgpm=list()
tm=list()

for(h in 1:5){
  
  funoc=function(x,i) {
    if     (i==1){5+2*cos(x)}
    else if(i==2){-5+2*sin(1*x)}
    
  }
  
  rl1=c()
  rl2=c()
  rl3=c()
  simulation=5000
  for(i in 1:simulation){
    for(j in 1:5000){
      ##Generating new Out-of-control profile ##
      yn=lapply(1:n,  function(i){funoc(trains,i)+rnorm(dp,0,nsd+h/100)#+c(rep(0,10),rep(h/10,10),rep(0,10))
        
      })
      t1=c()
      ##Forming T^2 Statistic##  Kriging  ###
      t1=unlist(lapply(1:n, function(i){mu1=0;mu1=yn[[i]]-ypredk[[i]];t(mu1)%*%solve(yvark[[i]],mu1)}))
      if(any(t1>q1)) {rl1[i]=j;break}
    }
    print(i)
  }
  
  for(i in 1:simulation){
    for(j in 1:5000){
      ##Generating new Out-of-control profile ##
      yn=lapply(1:n,  function(i){funoc(trains,i)+rnorm(dp,0,nsd+h/100)#+c(rep(0,10),rep(h/10,10),rep(0,10))
        
      })
      ##Forming T^2 Statistic## MGP ###
      mu2=unlist(yn)-ypred
      t2=as.numeric(t(mu2)%*%solve(yvar,mu2))
      if(t2>q2) {rl2[i]=j;break}
    }
    print(i)
  }
  
  for(i in 1:simulation){
    for(j in 1:5000){
      ##Generating new Out-of-control profile ##
      yn=lapply(1:n,  function(i){funoc(trains,i)+rnorm(dp,0,nsd+h/100)
      })
      t3=c()
      ##Forming T^2 Statistic##Typical T^2###
      t3=unlist(lapply(1:n, function(i){mu3=0;mu3=yn[[i]]-trainy[[i]];t(mu3)%*%solve(sy[[i]],mu3)}))
      if(any(t3>q3)) {rl3[i]=j;break}
    }
    print(i)
  }
  arl1=mean(rl1,na.rm=TRUE)
  arl2=mean(rl2,na.rm=TRUE)
  arl3=mean(rl3,na.rm=TRUE)
  
  gpm[h]=arl1
  mgpm[h]=arl2
  tm[h]=arl3
}
cat("Average Run Length of Kriging:")
print(gpm)
cat("Average Run Length of MGP:")
print(mgpm)
cat("Average Run Length of T^2:")
print(tm)
