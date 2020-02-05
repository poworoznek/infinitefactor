#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;

// [[Rcpp::export]]
Rcpp::NumericMatrix msf(arma::mat lambda, arma::mat pivot) {
  arma::mat refr = join_rows(lambda, -lambda);
  int k = lambda.n_cols;
  uvec ind = regspace<uvec> (0, k-1);
  uvec perm(k);
  arma::vec signs(k);
  rowvec norms(2*k);
  unsigned int w, c, wc;
  arma::mat diff, diffsq;
  
  for(int i=0; i<k; i++){
    diff = refr.each_col() - pivot.col(i);
    diffsq = square(diff);
    norms = sum(diffsq);
    w = index_min(norms);
    c = refr.n_cols / 2;
    if(w>=c){
      wc = w-c;
      signs(i) = -1;
      perm(i) = ind(wc);
      refr.shed_col(w);
      refr.shed_col(wc); 
      ind.shed_row(wc);} 
    else {
      wc = w+c;
      signs(i) = 1;
      perm(i) = ind(w);
      refr.shed_col(wc); 
      refr.shed_col(w);
      ind.shed_row(w);}
    }
  
  arma::mat permmat = zeros<arma::mat>(k,k);
  for(int i=0; i<k; i++){
    permmat(perm(i), i) = signs(i);
  }
  
  lambda *= permmat;
  return Rcpp::wrap(lambda);
}

// [[Rcpp::export]]
Rcpp::NumericVector msfOUT(arma::mat lambda, arma::mat pivot) {
  arma::mat refr = join_rows(lambda, -lambda);
  int k = lambda.n_cols;
  uvec ind = regspace<uvec> (0, k-1);
  uvec perm(k);
  arma::vec signs(k);
  rowvec norms(2*k);
  unsigned int w, c, wc;
  arma::mat diff, diffsq;
  
  for(int i=0; i<k; i++){
    diff = refr.each_col() - pivot.col(i);
    diffsq = square(diff);
    norms = sum(diffsq);
    w = index_min(norms);
    c = refr.n_cols / 2;
    if(w>=c){
      wc = w-c;
      signs(i) = -1;
      perm(i) = ind(wc);
      refr.shed_col(w);
      refr.shed_col(wc); 
      ind.shed_row(wc);} 
    else {
      wc = w+c;
      signs(i) = 1;
      perm(i) = ind(w);
      refr.shed_col(wc); 
      refr.shed_col(w);
      ind.shed_row(w);}
  }
  
  arma::vec out = (perm + ones<arma::vec>(k)) % signs;
  
  return Rcpp::wrap(out);
}

// [[Rcpp::export]]
Rcpp::NumericMatrix aplr(arma::mat matr, arma::vec perm){
  int k = matr.n_cols;
  arma::mat permmat = zeros<arma::mat>(k,k);
  arma::vec perms = abs(perm) - ones<arma::vec>(k);
  arma::vec signs = sign(perm);
  
  for(int i=0; i<k; i++){
    permmat(perms(i), i) = signs(i);
  }
  
  matr *= permmat;
  return(Rcpp::wrap(matr));
}
