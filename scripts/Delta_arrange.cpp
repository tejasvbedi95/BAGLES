// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp ;

static IntegerVector Vector_equal(IntegerVector x);
static IntegerVector Locator(IntegerVector x);

// [[Rcpp::export]]
IntegerVector delta_arrange( int T, int M ) {
  int t;
  int count = 0;
  int i;
  int sum = 0;
  int sum2 = 0;
  int k;
  IntegerVector delta_vec (T);
  for( t = 0; t < T; t+=ceil(T/M)+1 ){
    delta_vec(t) = 1;
    sum2 = sum2 + delta_vec(t);
    count++;
    if( sum2 > M ){
      delta_vec(t) = 0;
    }
  }
  IntegerVector delta_vec_temp = Vector_equal(delta_vec);
  for(i = T - 7; i < T; i++){
    if(delta_vec(i) == 1){
      delta_vec_temp(i) = 0;
      delta_vec_temp(T - 7) = 1;
    }
    sum = sum + delta_vec(i);
  }
  IntegerVector loc = Locator( delta_vec_temp );
  IntegerVector loc_diff = diff( loc );
  for (k = 0; k < loc_diff.length(); k++){
    if ( loc_diff[k] < 7){
      Rcerr << "Error : Too many change points or too few time points ";
      IntegerVector e1;
      return e1;
    }
  }
  if ( count < M || sum > 1){
    Rcerr << "Error : Too many change points or too few time points ";
    //std::exit(EXIT_FAILURE);
    IntegerVector e2;
    return e2;
  }
  return delta_vec_temp;
}

// [[Rcpp::export]]
IntegerVector Vector_equal(IntegerVector x) {
  int i;
  IntegerVector y (x.length());
  for(i = 0; i < x.length(); i++ ){
    y[i] = x[i];
  }
  return y;
}

// [[Rcpp::export]]
IntegerVector Locator(IntegerVector x) {
  IntegerVector loc (1);
  int t;
  for( t = 1; t < x.length(); t++ ){
    if( x[t] == 1 ){
      loc.push_back(t);
    }
  } 
  return loc;
}
