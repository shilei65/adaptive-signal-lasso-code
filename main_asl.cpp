//[[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
using namespace Rcpp;
using namespace Eigen; 
using namespace std;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

double rcpp_min(double x, double y) {
  double out;
  if (x < y) {
    out=x;
  }
  else{
    out=y;
  }
  return out;
}

double rcpp_max(double x, double y) {
  double out;
  if (x < y) {
    out=y;
  }
  else{
    out=x;
  }
  return out;
}

double rcpp_sign(double x){
  double out;
  if (x<0){
    out=-1;
  }
  else {
    out=1;
  }
  return out;
}

// [[Rcpp::export]]
List Signal_c(VectorXd Y, MatrixXd X, VectorXd beta0, VectorXd Re2, 
                  VectorXd ep1, VectorXd ep2,int constant, int iter_max, int p, double delta) {
  VectorXd beta1;
  VectorXd eps;
  double mu;
  int iter=1;
  double res=1, Re, Re1;
  while(iter<=iter_max && res>delta){
    beta1=beta0;
    eps=Y-X*beta0;
    if(constant==1){
      mu=eps.mean();
    }
    else{
      mu=0;
    }
    eps=eps.array()-mu;
    for(int i=0; i<p; i++){
       eps=eps+X.col(i)*beta0(i);
       Re1=eps.adjoint()*X.col(i);
       Re=Re1/Re2(i);
       if(Re<=0){
         beta0(i)=rcpp_min(0,(Re+ep1(i)));
       }
       if(Re>0 && Re<=1+ep2(i)){
         beta0(i)=rcpp_max(0,(Re-ep2(i)));
       }
       if(Re>1+ep2(i)){
         beta0(i)=rcpp_max(1,(Re-ep1(i))); 
       }
       eps=eps-X.col(i)*beta0(i);
    }
    iter=iter+1;
    res=sqrt((beta1-beta0).adjoint()*(beta1-beta0));
  }
  List out;
   out["Mu"]=mu;
   out["Beta"]=beta0;
   out["iters"]=iter;
  return out;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

