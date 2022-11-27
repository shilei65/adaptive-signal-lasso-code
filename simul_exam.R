
source("/Users/mac/Desktop/New\ signal\ lasso\ method/R\ code\ for\ adaptive\ SL/main_asl.R")
source("/Users/mac/Desktop/New\ signal\ lasso\ method/R\ code\ for\ adaptive\ SL/simulation_asl.R")
setwd("/Users/mac/Desktop/New\ signal\ lasso\ method/R\ code\ for\ adaptive\ SL")
Rcpp::sourceCpp("/Users/mac/Desktop/New\ signal\ lasso\ method/R\ code\ for\ adaptive\ SL/main_asl.cpp")

library(Rcpp)
library(lars)
library(msgps)
library(MASS)
library(Matrix)
library(glmnet)
library(ncvreg)

###caseI 1 2 3 gaussian####
n=100
p=30
p1=6
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(1)){
print(sigma)  
print(simulation_asl(n,beta,k=50,sigma,errortype="gaussian",disigntype="correlated",printout = FALSE))
}


n=50
p=150
p1=6
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(1,2)){
  print(sigma)  
  print(simulation_asl(n,beta,k=100,sigma,errortype="AR1",disigntype="correlated",printout = FALSE))
}


n=50
p=150
p1=20
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(1,2)){
  print(sigma)  
  print(simulation_asl(n,beta,k=100,sigma,errortype="AR1",disigntype="correlated",printout = FALSE))
}





###caseI 7 8 9 gaussian####
n=100
p=30
p1=6
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(1,2)){
  print(sigma) 
  print(simulation_asl(n,beta,k=100,sigma,errortype="gaussian",disigntype="uncorrelated",printout = FALSE))
}

###caseI 13 14 15 gaussian####
n=50
p=150
p1=6
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(1, 2)){
  print(sigma)  
  print(simulation_asl(n,beta,k=50,sigma,errortype="AR1",disigntype="correlated",printout = FALSE))
}


###caseI 19 20 21 gaussian####
n=50
p=150
p1=20
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(0.4,1)){
  print(sigma) 
  print(simulation_asl(n,beta,k=100,sigma,errortype="gaussian",disigntype="correlated",printout = FALSE))
} 

# ###caseI 22 23 24 gaussian####
n=30
p=200
p1=150
beta=c(rep(1,p1),rep(0,p-p1))
set.seed(1)
for(sigma in c(0.1)){
  print(sigma) 
  print(simulation_asl(n,beta,k=50,sigma,errortype="gaussian",disigntype="correlated",printout = FALSE))
} 
