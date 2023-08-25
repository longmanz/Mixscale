#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;


NumericVector vecpow(const NumericVector base, const IntegerVector exp) {
    NumericVector out(base.size());

    std::transform(base.cbegin(), base.cend(), exp.cbegin(), out.begin(), [](double a, double b) {return ::pow(a, b); });
    
    return out;
}


// [[Rcpp::export]]
double rbo_ext(StringVector x, StringVector y, double p, int k, bool uneven_lengths = true) {
    int l, s;
    StringVector S, L;
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
        NumericVector Xd(l);
        for(int i = 0; i < l; i++) {
            if(i < s){
               StringVector S_sub = S[Range(0, i)];  // Remember, indexing is 0-based in C++
               StringVector L_sub = L[Range(0, i)]; 
               Xd[i] = Rcpp::intersect(S_sub, L_sub).size();
            } else {
               StringVector L_sub = L[Range(0, i)]; 
               Xd[i] = Rcpp::intersect(S, L_sub).size();
            }
        }
        
        NumericVector rep1 = rep(p, l);
        NumericVector rep2 = rep(p, l-s);
        IntegerVector seq1 = seq_len(l);
        IntegerVector seq2 = seq_len(l-s);
        
        double sum1 = sum(Xd[seq(0, l-1)] / as<NumericVector>(seq1) * vecpow(rep1, seq1));
        double sum2 = sum(Xd[s-1] * (as<NumericVector>(seq2) + s) / (s * as<NumericVector>(seq2)) * vecpow(rep2, seq2));
        
        return ((1-p) / p) * (sum1 + sum2) + ((Xd[l-1] - Xd[s-1]) / l + (Xd[s-1] / s)) * pow(p, l);
        
    } else {
        //stopifnot(l == s);
        k = std::min(s, k);
        NumericVector Xd(k);
        for(int i = 0; i < k; i++) {
            StringVector S_sub = S[Range(0, i)];  // Remember, indexing is 0-based in C++
            StringVector L_sub = L[Range(0, i)]; 
            Xd[i] = Rcpp::intersect(S_sub, L_sub).size();
        }
        
        double Xk = Xd[k - 1];
        IntegerVector range = seq_len(k);
        
        NumericVector rep1 = rep(p, k);
        
        return (Xk / k) * pow(p, k) + ((1 - p) / p) * sum(Xd / as<NumericVector>(range) * vecpow(rep1, range));
    }
}




// [[Rcpp::export]]
double rbo_ext2(StringVector x, StringVector y, double p, int k, bool uneven_lengths = true) {
    int l, s;
    StringVector S, L;
    
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
    
    std::set<std::string> set_S(S.begin(), S.end()), set_L(L.begin(), L.end());
    
    NumericVector Xd(l);
    
    for(int i = 0; i < l; i++) {
        std::set<std::string> temp_S(set_S.begin(), std::next(set_S.begin(), std::min(i + 1, s)));
        std::set<std::string> temp_L(set_L.begin(), std::next(set_L.begin(), i + 1));
        
        std::vector<std::string> intersection;
        std::set_intersection(temp_S.begin(), temp_S.end(), temp_L.begin(), temp_L.end(), std::back_inserter(intersection));
        Xd[i] = intersection.size();
    }
    
    NumericVector rep1 = rep(p, l);
    NumericVector rep2 = rep(p, l-s);
    IntegerVector seq1 = seq_len(l);
    IntegerVector seq2 = seq_len(l-s);
    
    double sum1 = sum(Xd[seq(0, l-1)] / as<NumericVector>(seq1) * vecpow(rep1, seq1));
    double sum2 = sum(Xd[s-1] * (as<NumericVector>(seq2) + s) / (s * as<NumericVector>(seq2)) * vecpow(rep2, seq2));
    
    if (uneven_lengths) {
        return ((1-p) / p) * (sum1 + sum2) + ((Xd[l-1] - Xd[s-1]) / l + (Xd[s-1] / s)) * pow(p, l);
    } else {
        k = std::min(s, k);
        double Xk = Xd[k - 1];
        IntegerVector range = seq_len(k);
        
        NumericVector rep1 = rep(p, k);
        
        return (Xk / k) * pow(p, k) + ((1 - p) / p) * sum(Xd / as<NumericVector>(range) * vecpow(rep1, range));
    }
}




