#generate regression data

gen.regression.data=function(beta,X,n,p,sigma){
   y=c()
  
   y=  t(X)%*%beta +rnorm(n,mean=0,sd= sigma ) 
    
  
  
  return(y)
  
  
}
hausdorff.distance=function(v1,v2,n){
  p1=length(v1)
  p2=length(v2)
  if ( p1*p2==0){ dis=n} else{
    distance.mat=matrix(0,nrow=p1,ncol=p2)
    for( i in 1: p1){
      for( j in 1: p2){
        distance.mat[i,j]=abs(v1[i]-v2[j])
        
      }
    }
    dis=max(max( apply(distance.mat, 1, min)),max( apply(distance.mat, 2, min)))
  }
  return(dis)
  
}


colSD=function(matrix){
  Ncol=ncol(matrix)
  result=sapply(c(1:Ncol), function(t) sd(matrix[,t]))
  return(result)  
  
}





#random sample M intervals
gen.intervals=function(N,M){
  result.list=matrix(nrow=M,ncol=2)
  for( i in 1:M){
    result.list[i,]=c(sort(round(runif(2,min=0,max=N))))
    
  }
  return(result.list)
}