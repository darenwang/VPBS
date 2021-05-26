#cv for vpcusum

#compute the list of estimators
vpcusum.list= function(ran.intervals, y,X,delta ,lambda.group.list,vp.tau.list, y.new,X.new ){
  
   
  TT=length(vp.tau.list)
   LL= length(lambda.group.list)
   estimate.list=vector("list", length=LL)
   for (ll in 1: LL){
     print(ll)
     lambda.group=lambda.group.list[ll]
     record.temp=  cusum.record   ( ran.intervals  , y, X,delta ,lambda.group ,y.new,X.new) 
     estimate.list[[ll]]=    vector("list", length=TT)

     
  for (tt in 1: TT){
    est =  vpcusum ( N,ran.intervals, delta ,vp.tau.list[[tt]],record.temp    ) 
    if ( length(est)==0){
      estimate.list[[ll]][[tt]]=c(0,length(y))
    }else{
    
    estimate.list[[ll]][[tt]] = c(0, est,length(y))
    }
    } 
   }
   return(estimate.list)
}

#compute test errors from change points
test.error.from.change.point=function(lambda.group, y.train,X.train,y.test,X.test,estimation){
  res.total=0
  options(warn=-1)
  
  for ( kk in 2:length(estimation)){
    s=estimation[kk-1]
    e=estimation[kk]
    D=matrix(0, nrow=p,ncol=p)
    diag(D)=1
    
    out =  glmnet(t(X.train[,(s+1):e]), y.train[(s+1):e], family=c("gaussian"), 
                  alpha = 1,lambda=lambda.group/sqrt(e-s-1))
     res.temp= sum(( y.test[(s+1):e]- predict(out, newx=t(X.test[,(s+1):e]), lambda=lambda.group/sqrt(e-s+1)) )^2)
    #
     res.total=res.total+res.temp
  }
  
  return(res.total)
  
}

cv.vpcusum =function(  ran.intervals, y, X,delta ,lambda.group.list,vp.tau.list,y.new,X.new ){
  N=length(y)
 
   
  estimate.list=
    vpcusum.list( ran.intervals, y , X ,delta ,lambda.group.list,vp.tau.list,y.new,X.new )
  cv.matrix=matrix(0, nrow= length(lambda.group.list), ncol=length(vp.tau.list))
  
  for ( j in 1:length(lambda.group.list)){
    lambda.group=lambda.group.list[j]
    
    for ( ll in 1: length(vp.tau.list)){
      cv.matrix[j,ll]=test.error.from.change.point(lambda.group, y ,X ,y.new,X.new,
                                                   estimate.list[[j]][[ll]])
      
    }
  }
  index= which(cv.matrix == min(cv.matrix), arr.ind = TRUE)
  if(nrow(index)>1){index=index[1,]}
 # lambda.group.choose=lambda.group.list[index[1]]* 2 ^0.5
  #vp.tau.choose= vp.tau.list[index[2]]*2^0.5
  #new.intervals= gen.intervals(N,M)
 # record.choose= cusum.record   (new.intervals , y, X,delta ,lambda.group.choose)
 # result = vpcusum(new.intervals, y,X,delta ,vp.tau,record.choose   ) 
  
  
  
  # print(estimate.list[[index[1]]][[index[2]]])
  return(estimate.list[[index[1]]][[index[2]]])
  #return( 2* estimation.list[[index[1]]][[index[2]]] )
  
}
