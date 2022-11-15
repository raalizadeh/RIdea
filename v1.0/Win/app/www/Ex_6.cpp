#include <Rcpp.h>
  using namespace Rcpp;

// [[Rcpp::export]]
SEXP Ex_6(IntegerVector y, IntegerVector x){
IntegerVector UNI = union_(x,y);
IntegerVector INTER = intersect(x,y);
return Rcpp::List::create(UNI,INTER);
}
