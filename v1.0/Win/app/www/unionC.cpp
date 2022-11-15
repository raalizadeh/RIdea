#include <Rcpp.h>
#include <unordered_set>
#include <algorithm>
using namespace Rcpp;


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]

IntegerVector unionC(IntegerVector x, IntegerVector y) {
    int nx = x.size();
    int ny = y.size();

    IntegerVector tmp(nx + ny);

    std::sort(x.begin(),x.end()); // unique
    std::sort(y.begin(),y.end());

IntegerVector::iterator out_end = std::set_union(
x.begin(), x.end(), y.begin(), y.end(), tmp.begin()
);

int prev_value = 0;
IntegerVector out;
for (IntegerVector::iterator it = tmp.begin();
it != out_end; ++it){
if ((it != tmp.begin()) && (prev_value == *it)) continue;
out.push_back(*it);
prev_value = *it;
}
return out;
}
 