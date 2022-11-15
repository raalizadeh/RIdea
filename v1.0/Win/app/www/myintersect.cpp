#include <Rcpp.h>
using namespace Rcpp;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector myintersect(NumericVector x, NumericVector y) {
    std::vector<double> res;
    std::unordered_set<double> s(y.begin(), y.end());
    for (int i=0; i < x.size(); ++i) {
        auto f = s.find(x[i]);
        if (f != s.end()) {
            res.push_back(x[i]);
            s.erase(f);
        }
    }
    return Rcpp::wrap(res);
}