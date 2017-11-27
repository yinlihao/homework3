#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
arma::mat SVTcol(arma::mat X,double lambda){
  arma::mat U,V,Ud,Vd,Xhat;arma::vec s,sd;
  arma::svd(U,s,V,X);
  Ud=U.rows(find(s>0.000001));
  Vd=V.cols(find(s>0.000001));
  sd=s.elem(find(s>0.000001));
  for(int i=0;i<sd.n_elem;i++){
    if(sd(i)>lambda){sd(i)=sd(i)-lambda;}
    else{sd(i)=0;}
  }
  Xhat=Ud*diagmat(sd)*trans(Vd);
  return Xhat;
}
arma::mat SVTrow(arma::mat X,double lambda){
  arma::mat U,V,Ud,Vd,Xhat;arma::vec s,sd;
  arma::svd(U,s,V,X);
  Ud=U.cols(find(s>0.000001));
  Vd=V.rows(find(s>0.000001));
  sd=s.elem(find(s>0.000001));
  for(int i=0;i<sd.n_elem;i++){
    if(sd(i)>lambda){sd(i)=sd(i)-lambda;}
    else{sd(i)=0;}
  }
  Xhat=Ud*diagmat(sd)*trans(Vd);
  return Xhat;
}
arma::mat shedmat(arma::mat a,double p){
  int m;int n;
  m=a.n_rows;n=a.n_cols;
  mat b;b.randu(m,n);
  for(int i=0;i<m;i++){
    for(int j=0;j<n;j++){
      if(b(i,j)<=p){a(i,j)=0;}
    }
  }
  return a;
}
arma::mat projectmat(arma::mat a,arma::mat b){
  for(int i=0;i<b.n_rows;i++){
    for(int j=0;j<b.n_cols;j++){
      if(b(i,j)==0){a(i,j)=0;}
    }
  }
  return a;
}
// [[Rcpp::export]]
arma::mat APO(arma::mat X,double lambda,double p,double c){
  arma::mat Xhat,Xpart,Yhat,lastX;int m,n;m=X.n_rows;n=X.n_cols;double delta=1;
  Xpart=shedmat(X,p);
  Xhat.randu(m,n);
  do{
    Yhat=Xhat+delta*projectmat((Xpart-Xhat),Xpart);
    lastX=Xhat;
    if(m>n){Xhat=SVTrow(Yhat,lambda*delta);}
    else{Xhat=SVTcol(Yhat,lambda*delta);}
  } while (norm(lastX-Xhat)>c);
  return Xhat;
}





