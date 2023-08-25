#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector return1(NumericVector x, NumericVector y, int k) {
    // NumericVector Xd(k);
    // 
    // for(int i = 0; i < k; i++) {
    //     // Changed subsetting using seq(0, i) to a more direct subsetting 
    //     Xd[i] = Rcpp::intersect(x[Range(0, i)], y[Range(0, i)]).size();
    // }
    
    return Rcpp::intersect(x[Range(0, k)], y[Range(0, k)]);
}

