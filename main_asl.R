
####the correlation function######
Tpm=function(p,rho){
  A=matrix(0,p,p)
  for (i in 1:p) {
    for (j in 1:p) {
      if(i==j){A[i,j]=1}
      else{A[i,j]=rho^(abs(i-j))}
    }
    
  }
  return(A)
}
########The adaptiv signal lasso Function lambda1*\omega1|X|+lambda2*\omega2 |X-1|########
SIGNAL_l=function(Y,X,beta0,lambda1,lambda2,weights=1,constant=TRUE,iter_max,delta){
  #Y: the response vector
  #X: the covariate matrix
  #beta0: initial value of the estimated parameter
  #lambda1: the first penalty parameter
  #lambda2: the second penalty parameter
  #weights: if 1 is the signal lasso and a vector of length p is the adaptive signal lasso
  #iter_max: the maximal iterations
  #delta: control the accuracy 
  p=ncol(X)
  Re2=numeric(p)
  for (i in 1:p) {
    Re2[i]=t(X[,i])%*%X[,i]  
  } 
  lambda3=weights*lambda2
  ep1=(lambda1+lambda3)/Re2
  ep2=(lambda1-lambda3)/Re2
  ptm=proc.time()
  
  con=1*constant
  fit.c=Signal_c(Y=Y, X=X, beta0=beta0, Re2=Re2, 
                 ep1=ep1, ep2=ep2,constant = con, iter_max=iter_max, p=p, delta=delta)
  mse=mean((Y-fit.c$Mu-X%*%fit.c$Beta)^2)
  return(list(Mu=fit.c$Mu, Beta=fit.c$Beta,Mse=mse,Iter=fit.c$iters,Times=proc.time()-ptm))
}

#######The CV function for signal lasso####
CV.SIGNAL=function(Y,X,beta0,nfolds,nlambda1,alpha=seq(0.3,0.7,length.out=5),weights=1,constant=TRUE,iter_max,delta){
  n=length(Y)
  
  folds=cv.folds(n,nfolds)
  
  re=NULL
  lambda3=NULL
  for (lambda2 in nlambda1) {
    for (j in alpha) {
      se=0
      lambda1=j*lambda2
      for (k in 1:nfolds) {
        Y1=Y[as.vector(unlist(folds[-k]))]
        X1=X[as.vector(unlist(folds[-k])),]
        Y2=Y[as.vector(unlist(folds[k]))]
        X2=X[as.vector(unlist(folds[k])),]
        fit=SIGNAL_l(Y=Y1,X=X1,beta0=beta0,lambda1,lambda2,weights=weights,constant=constant,iter_max,delta)
        se=se+mean((Y2-fit$Mu-X2%*%fit$Beta)^2)
      }
      re=c(re,se)
      lambda3=rbind(lambda3,c(lambda1,lambda2))
    }
  }
  index=which.min(re)
  return(list(lambda.1se=lambda3[index,],Re=re))   
}
