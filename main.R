

rm(list = ls())
set.seed(1)


######package and library

library(glmnet)
library(gglasso)
 #setwd("/Users/darenw/Dropbox/vpcusum/code/scene1/")

source("functions.R")
 source("vpcusum.R")
source("cv-vpcusum.R")
 library("MASS")

###### tuning parameters
lambda.group.list=c(0.5,1,2   )
vp.tau.list=seq( 1, 50, 3)


#########
n=100
p=100
sp=10
sigma=1
M=50
set.seed(0)
beta.list=list()
snr.can= c( 1, 1.2,1.4, 1.6 )
RR=100
vp.haus= matrix(0,nrow=RR,ncol=length(snr.can))
vp.list=  vector("list", length = length(snr.can) )
  


for ( ss in 1:length(snr.can)){
  snr= snr.can[ss] 
  delta=10
  vp.list[[ss]]=vector("list", length=RR)
  
  
  for ( i in 1: RR){
    
    print(ss*1000+i)
     

     #beta is the population matrix
    beta.list=list()
    bet.tem=2*seq(1:p)%%2-1
    beta.list[[1]]=c(bet.tem[1:sp],rep(0,p-sp))*snr
    beta.list[[2]]=-1*beta.list[[1]]
    beta.list[[3]]=beta.list[[1]]
    #X is the design matrix
    cov_mat <- stats::toeplitz(0.6^(0:(p-1)))
     X1=matrix(mvrnorm(n, rep(0, p), Sigma = cov_mat),nrow=p)
    X2=matrix(mvrnorm(n, rep(0, p), Sigma = cov_mat),nrow=p)
    X3=matrix(mvrnorm(n, rep(0, p), Sigma = cov_mat),nrow=p)
     
    X.data=cbind(X1, X2, X3)
    y.data=c(gen.regression.data(beta.list[[1]],X1,n,p,sigma),
             gen.regression.data(beta.list[[2]],X2,n,p,sigma),
             gen.regression.data(beta.list[[3]],X3,n,p,sigma))
     
    N=length(y.data)
    #population quantities
     
    
    true.changes=c(0, n , 2*n,  N)
    K=length(beta.list)
    
    
    
    #VP cusum
    y.odd= y.data[seq(1,(N-1),2)]
    X.odd= X.data[,seq(1,(N-1),2)]
    y.even =  y.data[seq(2, N  ,2)]
    X.even= X.data[,seq(2, N  ,2)]
    ran.intervals=gen.intervals(N/2 ,M)
     vp.estimate = cv.vpcusum  (  ran.intervals, y.odd, X.odd,delta ,lambda.group.list,vp.tau.list,y.even,X.even )
    print("vp="); print( vp.estimate)
    
    vp.haus[ i,ss ] = hausdorff.distance(c(0,2*vp.estimate,N) ,true.changes,N)
    print(paste("vp.haus =", mean(vp.haus[1:i,ss])))
    vp.list[[ss]][[i]] = vp.estimate
   
    
    ###end of vp 
    
    
    
     
  }
}

 
colSD=function(matrix){
  Ncol=ncol(matrix)
  result=sapply(c(1:Ncol), function(t) sd(matrix[,t]))
  return(result)  
  
}

colMeans(vp.haus/N) 
colSD(vp.haus/N)
colMeans(bsa.haus/N)
colSD(bsa.haus/N)

colMeans(MDC.haus/N)
colSD(MDC.haus/N)

colMeans(old.haus)
colSD(old.haus)

mean(vp.time) 
sd(vp.time)
mean(bsa.time) 
sd(bsa.time)
mean(sgl.time)
sd(sgl.time)
