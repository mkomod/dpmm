#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat rdirchlet(int n, arma::vec a)
{
    int k = a.n_elem;
    arma::mat res = arma::mat(n, k, arma::fill::zeros);
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < k; ++j) {
	    res(i, j) = R::rgamma(a(j), 1);
	}
	res.row(i) /= sum(res.row(i));
    }
    return res;
}

// [[Rcpp::export]]
arma::mat rdiscrete(int n, arma::vec p, bool as_indicator=false)
{
    int k = p.n_elem;
    arma::vec p_cs = arma::cumsum(p);
    arma::mat res;

    if (as_indicator) {
	res = arma::mat(n, k, arma::fill::zeros);
	for (int i = 0; i < n; ++i) {
	    int j = arma::find(p_cs > R::runif(0, 1), 1).eval()(0);
	    res(i, j) = 1;
	}
    } else {
	res = arma::mat(n, 1, arma::fill::zeros);
	for (int i = 0; i < n; ++i) {
	    int j = arma::find(p_cs > R::runif(0, 1), 1).eval()(0);
	    res(i) = j;
	}
    }
    return res;
}
