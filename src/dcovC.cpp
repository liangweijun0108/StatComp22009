#include <Rcpp.h>
using namespace Rcpp;

//' @title Computing brownian distance covariance using Rcpp
//' @description Computing brownian distance covariance of the input samples using Rcpp
//' @param x a group of samples
//' @param y a group of samples
//' @return the brownian distance covariance of the two groups of samples
//' @examples
//' \dontrun{
//' data(Eckerle4)
//' attach(Eckerle4)
//' dcov <- dcovC(Eckerle4$wavelength, Eckerle4$transmitance)
//' }
//' @export
// [[Rcpp::export]]
double dcovC(NumericVector x, NumericVector y) {
  int n = x.size();
  NumericMatrix a(n, n), b(n, n), A(n, n), B(n, n);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      a(i,j) = abs(x[i] - x[j]);
      b(i,j) = abs(y[i] - y[j]);
    }
  }
  NumericVector a1 = rowMeans(a);
  NumericVector a2 = colMeans(a);
  NumericVector b1 = rowMeans(b);
  NumericVector b2 = colMeans(b);
  for(int i = 0; i < n; i++){
    for(int j = 0; j < n; j++){
      A(i,j) = a(i,j) + mean(a) - a1[i] - a2[j];
      B(i,j) = b(i,j) + mean(b) - b1[i] - b2[j];
    }
  }
  return (sqrt(mean(A*B)));
}
