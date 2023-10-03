#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
arma::mat pw_l1_norm(arma::mat v1) {
  
  arma::mat l1_norm = arma::zeros<arma::mat>(v1.n_rows, v1.n_rows);
  
  for(int i=0;i< v1.n_rows; i++){
    for(int j = 0; j <= i; j++){
      l1_norm(i,j) = 1.0/float(v1.n_cols) * arma::sum(arma::abs( v1.row(i) - v1.row(j)));
    }
  }
  l1_norm = arma::symmatl(l1_norm);
  return(1.0 - l1_norm);
}


// [[Rcpp::export]]
float NDC(arma::mat v1, arma::mat v2){
  arma::mat e1 = pw_l1_norm(v1); 
  arma::mat e2 = pw_l1_norm(v2); 
  
  arma::uvec inds = arma::trimatl_ind(arma::size(e1));
  
  arma::vec e1v = e1(inds);
  arma::vec e2v = e2(inds);
  
  float ndc = 1.0 - arma::mean(arma::abs(e1v-e2v));
  return(ndc);
}

// [[Rcpp::export]]
float ENDC(arma::mat v1, arma::mat v2, int iter = 1000){
  arma::vec temp = arma::zeros<arma::vec>(iter);
  arma::mat e1 = pw_l1_norm(v1); 
  arma::mat e2 = pw_l1_norm(v2); 
  arma::uvec inds = arma::trimatl_ind(arma::size(e1));
  arma::vec e1v = e1(inds);
  arma::vec e2v = e2(inds);
  
  for(int i = 0; i < iter; i++){
   arma::vec v1s = arma::shuffle(e1v);
   arma::vec v2s = arma::shuffle(e2v);
   temp(i) = 1.0 - arma::mean(arma::abs(v1s-v2s));
  }
  return(float(arma::mean(temp)));
}

// [[Rcpp::export]]
float ACI(arma::mat v1, arma::mat v2, int iter = 1000){
  float ndc = NDC(v1, v2);
  float endc = ENDC(v1, v2, iter = iter);
  
  return((ndc-endc)/(1.0-endc));
  
}

// [[Rcpp::export]]
arma::mat rand_corr(int d, arma::mat rands) {
  
  
  arma::mat rand_corr = arma::eye(d, d);
  rand_corr.diag(1) = rands.diag(1);
  rand_corr.diag(-1) = rands.diag(1);
  for(int k = 2; k < (d); k++){
    for(int j = 0; j < (d-k); j++){
      // Rcout << j+1 << " " << j+k+1 <<std::endl;
      arma::rowvec r_1 = rand_corr(j, arma::span(j+1, j+k-1));
      //  Rcout << r_1 << std::endl;
      arma::rowvec r_3 = rand_corr(j+k, arma::span(j+1, j+k-1));
      //  Rcout << r_3 << std::endl;
      arma::mat R_2 = rand_corr.submat(j+1, j+1, j+k-1, j+k-1);
      //  Rcout << R_2 << std::endl;
      double rand = rands(j,j+k);
      
      arma::mat D_2jk = (1 - r_1*arma::inv(R_2)*arma::trans(r_1))*(1 - r_3*arma::inv(R_2)*arma::trans(r_3));
      rand_corr(j,j+k) = arma::as_scalar(r_1*arma::inv(R_2)*arma::trans(r_3) + rand*D_2jk);
      rand_corr(j+k,j) =rand_corr(j,j+k);
    }
  }
  return(rand_corr);
}

