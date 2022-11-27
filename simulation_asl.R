AR.1=function(n,rho,sigma){
  x=numeric(n)
  x[1]=rnorm(1,sd=sigma)
  for(t in 2:n){
    x[t]=rho*x[t-1]+rnorm(1,sd=sigma)
  }
  return(x)
}

##########################


error_compute<-function(beta, beta1){
  signal.1=1*I(beta1>0.9 & beta1<1.1)
  signal.0=1*I(beta1>-0.1 & beta1<0.1)
  TP=sum(signal.1*beta)
  FP=sum(signal.1*(1-beta))
  TN=sum(signal.0*(1-beta))
  FN=sum(signal.0*beta)
  UC=length(beta)-TP-FN-FP-TN
  UC1=sum(beta)-TP-FN
  UC0=sum(1-beta)-TN-FP
  return(list(TP = TP, FP = FP, TN = TN, FN = FN, UC = UC, UC1 = UC1, UC0 = UC0))
}


#######simulation function####
simulation_asl=function(n,beta,k,sigma,errortype="gaussian",disigntype="uncorrelated",printout=TRUE){
#n: number of samples
#beta: parameters  
#k: repeat times 
#sigma: scale parameters for different distributions
#errortype: include "gaussian", "exponential","gamma","AR(1)"
#disigntype: include "uncorrelated","correlated" 
#printout: print the iter times
p=length(beta)
p1=sum(I(beta!=0))

beta.l=TPR.l=TNR.l=MCC.l=FPR.l=SREL.l=SRNL.l=PCR.l=MCCs.l=G.l=UC.l=NULL  ##lasso
beta.al=TPR.al=TNR.al=MCC.al=FPR.al=SREL.al=SRNL.al=PCR.al=MCCs.al=G.al=UC.al=NULL  ##adaptive lasso
beta.sc=TPR.sc=TNR.sc=MCC.sc=FPR.sc=SREL.sc=SRNL.sc=PCR.sc=MCCs.sc=G.sc=UC.sc=NULL  ##SCAD lasso
beta.mc=TPR.mc=TNR.mc=MCC.mc=FPR.mc=SREL.mc=SRNL.mc=PCR.mc=MCCs.mc=G.mc=UC.mc=NULL  ##MCP lasso
beta.el=TPR.el=TNR.el=MCC.el=FPR.el=SREL.el=SRNL.el=PCR.el=MCCs.el=G.el=UC.el=NULL  ##elastic net
beta.s=TPR.s=TNR.s=MCC.s=FPR.s=SREL.s=SRNL.s=PCR.s=MCCs.s=G.s=UC.s=NULL  ##signal lasso
beta.as=TPR.as=TNR.as=MCC.as=FPR.as=SREL.as=SRNL.as=PCR.as=MCCs.as=G.as=UC.as=NULL ##adaptive signal lasso


for(t in 1:k){
  if(disigntype=="uncorrelated"){
  X=matrix(rnorm(n*p),n,p)
  }
  if(disigntype=="correlated"){
    X=mvrnorm(n,mu=rep(0,p),Sigma = Tpm(p,0.5))
  }
  if(errortype=="gaussian"){
  Y=X%*%beta+rnorm(n,sd=sigma)
  }
  if(errortype=="exponential"){
    Y=X%*%beta+rexp(n,rate=1/sigma)
  }
  if(errortype=="gamma"){
    Y=X%*%beta+rgamma(n,shape = 4,scale=sigma/2)
  }
  if(errortype=="AR1"){
    Y=X%*%beta+AR.1(n,rho=0.8,sigma)
  }
  

  ###lasso fit#####
  cvfit0=cv.glmnet(x=X,y=Y,nfolds=5,intercept=TRUE,type.measure = "mse")
  # fit.l=glmnet(x=X,y=Y,lambda =cvfit0$lambda.1se,alpha = 1,intercept = TRUE)
  # beta1=coef.glmnet(fit.l)[-1]
  fit.l=msgps(X=X,y=as.vector(Y), penalty="enet", alpha=0, lambda=cvfit0$lambda.1se, intercept=TRUE)
  beta1=coef(fit.l)[,3][-1]
  
  beta.l=rbind(beta.l, beta1)
  
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.l<-c(SREL.l, TP/sum(beta))
  SRNL.l<-c(SRNL.l, TN/sum(1-beta))         
  TPR.l<-c(TPR.l, TP/(TP+FN+0.001))
  FPR.l<-c(FPR.l, FP/(TN+FP+0.001))
  MCC.l<-c(MCC.l, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.l<-c(PCR.l, TP/(TP+FP))
  # F.l<-c(F.l, 2*((TP/(TP+FP))*(TP/(TP+FN)))/((TP/(TP+FP))+(TP/(TP+FN))))
  # G.l<-c(G.l, sqrt((TP/(TP+FP))*(TP/(TP+FN))))
  UC.l<-c(UC.l, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.l<-c(MCCs.l, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  ####adaptive lasso####
  
  fit.al=msgps(X=X,y=as.vector(Y), penalty="alasso", gamma=1, lambda=cvfit0$lambda.1se, intercept=TRUE)
  beta1=coef(fit.al)[,3][-1]
  beta.al=rbind(beta.al, beta1)
  beta0=beta1
 
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.al<-c(SREL.al, TP/sum(beta))
  SRNL.al<-c(SRNL.al, TN/sum(1-beta))         
  TPR.al<-c(TPR.al, TP/(TP+FN+0.001))
  FPR.al<-c(FPR.al, FP/(TN+FP+0.001))
  MCC.al<-c(MCC.al, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.al<-c(PCR.al, TP/(TP+FP))
  # F.al<-c(F.al, 2*((TP/(TP+FP))*(TP/(TP+FN)))/((TP/(TP+FP))+(TP/(TP+FN))))
  # G.al<-c(G.al, sqrt((TP/(TP+FP))*(TP/(TP+FN))))
  UC.al<-c(UC.al, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.al<-c(MCCs.al, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  
  #### SCAD ####
  cvfit=cv.ncvreg(X=X,y=Y,nfolds=5,penalty="SCAD",intercept = TRUE)
  fit.sc=cvfit$fit
  beta1=fit.sc$beta[,cvfit$min][-1]
  beta.sc=rbind(beta.sc, beta1)
  
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.sc<-c(SREL.sc, TP/sum(beta))
  SRNL.sc<-c(SRNL.sc, TN/sum(1-beta))         
  TPR.sc=c(TPR.sc,TP/(TP+FN+0.001))
  FPR.sc<-c(FPR.sc, FP/(TN+FP+0.001))
  MCC.sc<-c(MCC.sc, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.sc<-c(PCR.sc, TP/(TP+FP))

  UC.sc<-c(UC.sc, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.sc<-c(MCCs.sc, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  #### MCP ####
  cvfit=cv.ncvreg(X=X,y=Y,nfolds=5,penalty="MCP",intercept = TRUE)
  fit.mc=cvfit$fit
  beta1=fit.mc$beta[,cvfit$min][-1]
  beta.mc=rbind(beta.mc, beta1)
  
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.mc<-c(SREL.mc, TP/sum(beta))
  SRNL.mc<-c(SRNL.mc, TN/sum(1-beta))         
  TPR.mc=c(TPR.mc,TP/(TP+FN))
  FPR.mc<-c(FPR.mc, FP/(TN+FP))
  MCC.mc<-c(MCC.mc, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.mc<-c(PCR.mc, TP/(TP+FP))
  
  UC.mc<-c(UC.mc, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.mc<-c(MCCs.mc, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  #### elasticnet ####
  lam=cvm=NULL
  alpha0=seq(0.1,0.9,0.1)
  for(alpha1 in seq(0.1,0.9,0.1)){
  cvfit=cv.glmnet(x=X,y=Y,nfolds=5,intercept=TRUE,type.measure = "mse",alpha=alpha1)
  lam=c(lam,cvfit$lambda.1se)
  cvm=c(cvm, min(cvfit$cvm))
  }
  fit.el=glmnet(x=X,y=Y,lambda =lam[which.min(cvm)],alpha =alpha0[which.min(cvm)],intercept = TRUE)
  beta1=coef(fit.el)[-1]
  beta.el=rbind(beta.el, beta1)
  
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.el<-c(SREL.el, TP/sum(beta))
  SRNL.el<-c(SRNL.el, TN/sum(1-beta))         
  TPR.el<-c(TPR.el, TP/(TP+FN))
  FPR.el<-c(FPR.el, FP/(TN+FP))
  MCC.el<-c(MCC.el, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.el<-c(PCR.el, TP/(TP+FP))
  
  
  UC.el<-c(UC.el, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.el<-c(MCCs.el, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  #### signal lasso####
  #fit=CV.SIGNAL(Y=Y,X=X,beta0 = beta0,nfolds=5,nlambda1=c(1,0.1,0.01,0.001,0.0001),alpha=c(seq(0.1,0.5,0.05)), 
  #              weights = 1,constant=TRUE,iter_max = 2000,delta=1e-7)
  
  fit=CV.SIGNAL(Y=Y,X=X,beta0 = beta0,nfolds=5,nlambda1=cvfit0$lambda[1:10],alpha=c(seq(0.1,0.9,0.1),1/seq(0.1,0.9,0.1)), 
                weights = 1,constant=TRUE,iter_max = 2000,delta=1e-7)
  fit.s=SIGNAL_l(Y=Y,X=X,beta0 = beta0,lambda1=fit$lambda.1se[1],
                 lambda2 =fit$lambda.1se[2],weights=1,iter_max = 2000,delta=1e-7)
  
  beta1=fit.s$Beta
  beta.s=rbind(beta.s, beta1)
  
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.s<-c(SREL.s, TP/sum(beta))
  SRNL.s<-c(SRNL.s, TN/sum(1-beta))         
  TPR.s<-c(TPR.s, TP/(TP+FN+0.0001))
  FPR.s<-c(FPR.s, FP/(TN+FP+0.0001))
  MCC.s<-c(MCC.s, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.s<-c(PCR.s, TP/(TP+FP))
  # F.s<-c(F.s, 2*((TP/(TP+FP))*(TP/(TP+FN)))/((TP/(TP+FP))+(TP/(TP+FN))))
  # G.s<-c(G.s, sqrt((TP/(TP+FP))*(TP/(TP+FN))))
  UC.s<-c(UC.s, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.s<-c(MCCs.s, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  
  ####adaptive signal lasso####
  fit=glmnet(x=X,y=Y,lambda =cvfit0$lambda.1se,alpha =0,intercept = TRUE)
  w=abs(coef.glmnet(fit)[-1])
  fit=CV.SIGNAL(Y=Y,X=X,beta0 = beta0,nfolds=5,nlambda1=2000,alpha=c(seq(0,0.8,0.01)), 
                weights =w,constant=TRUE,iter_max = 2000,delta=1e-7)
  fit.as=SIGNAL_l(Y=Y,X=X,beta0 = beta0,lambda1=fit$lambda.1se[1],
                 lambda2 =fit$lambda.1se[2],weights=w,iter_max = 2000,delta=1e-7)
  beta1=fit.as$Beta
  beta.as=rbind(beta.as, beta1)
  # 
  eps<-error_compute(beta, beta1)
  TP<-eps[["TP"]]; FP<-eps[["FP"]]; TN<-eps[["TN"]]; FN<-eps[["FN"]]; UC<-eps[["UC"]]; UC1<-eps[["UC1"]]; UC0<-eps[["UC0"]]
  SREL.as<-c(SREL.as, TP/sum(beta))
  SRNL.as<-c(SRNL.as, TN/sum(1-beta))         
  TPR.as<-c(TPR.as, TP/(TP+FN+0.001))
  FPR.as<-c(FPR.as, FP/(TN+FP+0.001))
  MCC.as<-c(MCC.as, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  PCR.as<-c(PCR.as, TP/(TP+FP))
  # F.as<-c(F.as, 2*((TP/(TP+FP))*(TP/(TP+FN)))/((TP/(TP+FP))+(TP/(TP+FN))))
  # G.as<-c(G.as, sqrt((TP/(TP+FP))*(TP/(TP+FN))))
  UC.as<-c(UC.as, UC/length(beta))
  FN<-FN+UC1
  FP<-FP+UC0
  MCCs.as<-c(MCCs.as, (TP*TN-FP*FN)/(sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))))
  
  
  if(printout){
  print(t)
    }
}
MCC.l[which(is.na(MCC.l))]=0
MCC.al[which(is.na(MCC.al))]=0
MCC.sc[which(is.na(MCC.sc))]=0
MCC.mc[which(is.na(MCC.mc))]=0
MCC.el[which(is.na(MCC.el))]=0
MCC.s[which(is.na(MCC.s))]=0
MCC.as[which(is.na(MCC.as))]=0


PCR.l[which(is.na(PCR.l))]=0
PCR.al[which(is.na(PCR.al))]=0
PCR.sc[which(is.na(PCR.sc))]=0
PCR.mc[which(is.na(PCR.mc))]=0
PCR.el[which(is.na(PCR.el))]=0
PCR.s[which(is.na(PCR.s))]=0
PCR.as[which(is.na(PCR.as))]=0


SREL.l[which(is.na(SREL.l))]=0
SREL.al[which(is.na(SREL.al))]=0
SREL.sc[which(is.na(SREL.sc))]=0
SREL.mc[which(is.na(SREL.mc))]=0
SREL.el[which(is.na(SREL.el))]=0
SREL.s[which(is.na(SREL.s))]=0
SREL.as[which(is.na(SREL.as))]=0


TPR.l[which(is.na(TPR.l))]=0
TPR.al[which(is.na(TPR.al))]=0
TPR.sc[which(is.na(TPR.sc))]=0
TPR.mc[which(is.na(TPR.mc))]=0
TPR.el[which(is.na(TPR.el))]=0
TPR.s[which(is.na(TPR.s))]=0
TPR.as[which(is.na(TPR.as))]=0


FPR.l[which(is.na(FPR.l))]=0
FPR.al[which(is.na(FPR.al))]=0
FPR.sc[which(is.na(FPR.sc))]=0
FPR.mc[which(is.na(FPR.mc))]=0
FPR.el[which(is.na(FPR.el))]=0
FPR.s[which(is.na(FPR.s))]=0
FPR.as[which(is.na(FPR.as))]=0

MCCs.l[which(is.na(MCCs.l))]=0
MCCs.al[which(is.na(MCCs.al))]=0
MCCs.sc[which(is.na(MCCs.sc))]=0
MCCs.mc[which(is.na(MCCs.mc))]=0
MCCs.el[which(is.na(MCCs.el))]=0
MCCs.s[which(is.na(MCCs.s))]=0
MCCs.as[which(is.na(MCCs.as))]=0


#####results####
Mse=c(mean((colMeans(beta.l)-beta)^2)+mean(apply(beta.l, 2, var)),
mean((colMeans(beta.al)-beta)^2)+mean(apply(beta.al, 2, var)),
mean((colMeans(beta.sc)-beta)^2)+mean(apply(beta.sc, 2, var)),
mean((colMeans(beta.mc)-beta)^2)+mean(apply(beta.mc, 2, var)),
mean((colMeans(beta.el)-beta)^2)+mean(apply(beta.el, 2, var)),
mean((colMeans(beta.s)-beta)^2)+mean(apply(beta.s, 2, var)),
mean((colMeans(beta.as)-beta)^2)+mean(apply(beta.as, 2, var)))


TPR=c(mean(TPR.l),mean(TPR.al),mean(TPR.sc),mean(TPR.mc),
      mean(TPR.el),mean(TPR.s),mean(TPR.as))
MCC=c(mean(MCC.l),mean(MCC.al),mean(MCC.sc),mean(MCC.mc),
      mean(MCC.el),mean(MCC.s),mean(MCC.as))
PCR=c(mean(PCR.l),mean(PCR.al),mean(PCR.sc),mean(PCR.mc),
      mean(PCR.el),mean(PCR.s),mean(PCR.as))


SREL=c(mean(SREL.l),mean(SREL.al),mean(SREL.sc),mean(SREL.mc),
      mean(SREL.el),mean(SREL.s),mean(SREL.as))
SRNL=c(mean(SRNL.l),mean(SRNL.al),mean(SRNL.sc),mean(SRNL.mc),
       mean(SRNL.el),mean(SRNL.s),mean(SRNL.as))
FPR=c(mean(FPR.l),mean(FPR.al),mean(FPR.sc),mean(FPR.mc),
      mean(FPR.el),mean(FPR.s),mean(FPR.as))

UCR=c(mean(UC.l),mean(UC.al),mean(UC.sc),mean(UC.mc),
      mean(UC.el),mean(UC.s),mean(UC.as))

MCCs=c(mean(MCCs.l),mean(MCCs.al),mean(MCCs.sc),mean(MCCs.mc),
       mean(MCCs.el),mean(MCCs.s),mean(MCCs.as))

out=rbind(Mse,SREL,SRNL,TPR,PCR,FPR,MCC,MCCs,UCR)
colnames(out)=c("lasso","a-lasso","SCAD","MCP",
                "elastic","signal","a-signal")
return(t(out))
}


