#one change point

#group lasso
convert.design.matrix.one.change=function(X, s,e,t){
  xx1=t(X[,(s+1):e])
  xx1[(t-s+1):(e-s),]=0
  
  xx2=t(X[,(s+1):e])
  xx2[ 1 :(t-s),]=0
  

  xx=cbind(xx1/sqrt(t-s),xx2/sqrt(e-t))
  
  
  return(xx)
}

one.change.group.lasso =function(t,s,e,y, X ,lambda.group){
  group=c(seq(1:p),c(1:p))
  out=gglasso(x=convert.design.matrix.one.change(X,s,e,t),y=y[(s+1):e],group=group, loss="ls",
              lambda=lambda.group/ (e-s),intercept = FALSE,eps = 0.0001  )
  
 
  alpha=as.vector(out$beta)
   beta=c(alpha[1:p ]/sqrt(t-s),alpha[(p+1):(2*p)]/sqrt(e-t))
  res=sum((y[(s+1):e]-convert.design.matrix.one.change(X,s,e,t) %*%alpha)^2)
  return( list("res"=res,"beta"=beta))
}
find.one.change.grouplasso=function(y,X,s,e,delta ,p ,lambda.group){
  can.vec=c((s+delta): (e-delta))
  res.seq=sapply(can.vec,function(t) 
    one.change.group.lasso(t,s,e,y, X ,lambda.group )$res)
  return(which.min(res.seq)+s+delta-1)
  
}

########cusum 1D
cusum.factor= function (t,s,e){
  return (( e-t)*(t- s) / (e-s) )
  
}
cusum.one.dimension=function(time.series,s,e,t,delta){
 cusum.temp=   (cusum.factor(t,s+1,e) )^0.5 * 
           ( mean(time.series[(s+1):t])  -mean(time.series[(t+1):e]) )   
return( abs(  cusum.temp ))
}

  


#####vpcusum
vpcusum.one.interval=function(s,e, y, X,delta ,lambda.group, y.new, X.new ){
  NN=length(y)
  estimate=find.one.change.grouplasso (y,X,s,e,delta ,p ,lambda.group)
  beta.temp = one.change.group.lasso  (estimate,s,e,y, X ,lambda.group)$beta
  projection.vector=     beta.temp[1:p] -beta.temp[(p+1):(2*p)] 
   if(sum(abs(projection.vector))==0){
      projection.vector = projection.vector  }else{
   projection.vector= projection.vector/norm(projection.vector, type='2')}
  projected.data=sapply (  1  : NN, function(x) y.new [x] *projection.vector%*%X.new [,x]   )
  return(list( "value"=cusum.one.dimension (projected.data,s,e,estimate), "estimate"=estimate))
}

#compute the record
 cusum.record =function ( ran.intervals  , y, X,delta ,lambda.group,y.new, X.new ){
  M=nrow(ran.intervals)
  NN=length(y)
  record=matrix(0,nrow=M, ncol=2)
  for (m in 1:M){
    #print(m)
    s=ran.intervals[m,1]
    e=ran.intervals[m,2]
    
     if(e -s >2*delta){
         out.one.interval =vpcusum.one.interval(s,e, y, X,delta ,lambda.group, y.new, X.new)
         record[m,1]= out.one.interval$value 
         record[m,2]= out.one.interval$estimate
    
    
     } 
  }
   return(record )
}

 
##########use the cusum record to compute change point
 
vpcusum=function( NN,ran.intervals, delta ,vp.tau,record    ){
 
  start.set=c(0,NN)
  output.set=c()
  while(length(start.set)>1){
    s=start.set[1]
    e=start.set[2]
    # print(start.set)
    
    
    record.temp= matrix(record[which(ran.intervals[ , 1]>s & ran.intervals[ ,2]<e+1) ,],ncol=2)
    #print(record.temp)
    if(nrow(record.temp)==0){start.set=start.set[-1]}else{
      out.temp= record.temp[which.max(record.temp [,1]),] 
      
      
      if(out.temp[ 1]>vp.tau){
        start.set=sort(c(start.set,out.temp[ 2]))
        output.set = (c(output.set,out.temp[ 2]))
      }else{
        start.set=start.set[-1]
      }
      
    }
  }
   return( sort(output.set))
  
  
}

