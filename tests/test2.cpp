#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector vecpow(const NumericVector base, const IntegerVector exp) {
    NumericVector out(base.size());
    
    std::transform(base.cbegin(), base.cend(), exp.cbegin(), out.begin(), [](double a, double b) {return ::pow(a, b); });
    
    return out;
}


// [[Rcpp::export]]
double rbo_ext(NumericVector x, NumericVector y, double p, int k, bool uneven_lengths = true) {
    int l, s;
    NumericVector S, L;
    if (x.length() <= y.length()) {
        S = x;
        L = y;
    } else {
        S = y;
        L = x;
    }
    
    // l = min(k, L.length());
    if(k <= L.length()){
        l = k;
    } else {
        l = L.length();
    }
    
    // s = min(k, S.length());
    if(k <= S.length()){
        s = k;
    } else {
        s = S.length();
    }
    
    if (uneven_lengths) {
       // NumericVector Xd(l);
       // for(int i = 0; i < l; i++) {
       //     if(i < s){
       //        NumericVector S_sub = S[Range(0, i)];  // Remember, indexing is 0-based in C++
       //        NumericVector L_sub = L[Range(0, i)]; 
       //        Xd[i] = Rcpp::intersect(S_sub, L_sub).size();
       //     } else {
       //        NumericVector L_sub = L[Range(0, i)]; 
       //        Xd[i] = Rcpp::intersect(S, L_sub).size();
       //     }
       // }
        
        // IntegerVector range_1 = seq_len(l);
        // IntegerVector range_2 = seq_len(l) - s;
        // NumericVector rep1 = rep(p, range_1.length());
        // NumericVector rep2 = rep(p, range_2.length());
        //     
        // NumericVector p_seq_1 = vecpow(rep1, range_1);
        // NumericVector p_seq_2 = vecpow(rep2, range_2 + s + 1);
        // 
        // double sum1 = sum(Xd / as<NumericVector>(range_1) * p_seq_1);
        // double sum2 = sum(Xd[s - 1] * as<NumericVector>(range_2) / (s * as<NumericVector>(range_2 + s + 1)) * p_seq_2);
        // 
        // return ((1 - p) / p) * (sum1 + sum2) + ((Xd[l - 1] - Xd[s - 1]) / l + (Xd[s - 1] / s)) * pow(p, l);
        
        
        //NumericVector rep1 = rep(p, l);
        //NumericVector rep2 = rep(p, l-s);
        //IntegerVector seq1 = seq_len(l);
        //IntegerVector seq2 = seq_len(l-s);
        //
        //double sum1 = sum(Xd[seq(0, l-1)] / as<NumericVector>(seq1) * vecpow(rep1, seq1));
        //double sum2 = sum(Xd[s-1] * (as<NumericVector>(seq2) + s) / (s * as<NumericVector>(seq2)) * vecpow(rep2, seq2));
        //
        //double result = ((1-p) / p) * (sum1 + sum2) + ((Xd[l-1] - Xd[s-1]) / l + (Xd[s-1] / s)) * pow(p, l);
        //return result;
        
        return 1.0;
        
    } else {
        //stopifnot(l == s);
        k = std::min(s, k);
        NumericVector Xd(k);
        for(int i = 0; i < k; i++) {
           NumericVector S_sub = S[Range(0, i)];  // Remember, indexing is 0-based in C++
           NumericVector L_sub = L[Range(0, i)]; 
           Xd[i] = Rcpp::intersect(S_sub, L_sub).size();
        }
        // double Xk = Xd[k - 1];
        // IntegerVector range = seq_len(k);
        
        // NumericVector rep1 = rep(p, k);
        
        // double result = (Xk / k) * pow(p, k) + ((1 - p) / p) * sum(Xd / as<NumericVector>(range) * vecpow(rep1, range));
        // return result;
        
        return 1.0; 
    }
}



