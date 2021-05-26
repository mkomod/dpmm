#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using arma::uword;

// [[Rcpp::export]]
double log_marginal_lik(double y, double x, double eta, double sigma, double sigma_0)
{
    double z = y - eta;
    double s = 2.0*sigma*sigma;
    double s_0 = 1.0/(2.0*sigma_0*sigma_0);
    double A = z * x / s;
    double B = x * x / s;
    double D = z * z / s;

    return (-D + A*A/(B + s_0));
}


double marginal_lik(double y, double x, double eta, double sigma, double sigma_0)
{
    return exp(log_marginal_lik(y, x, eta, sigma, sigma_0));
}


double marginal_lik_ratio(double eta_new, double eta_old, double y, double x,
	double sigma, double sigma_0)
{

    double res = marginal_lik(y, x, eta_new, sigma, sigma_0) /
		 marginal_lik(y, x, eta_old, sigma, sigma_0);
    return res;

}

// checks if the jth element in a is a singleton (i.e. unique)
bool is_singleton(const arma::uvec &a, uword j) 
{
    for (uword i = 0; i < a.n_elem; ++i) {
	if (i == j) continue;
	if (a(i) == a(j)) return false;
    }
    return true;
}


int num_groups(const arma::vec eta) 
{
    int n = eta.n_elem;
    for (int i = 0; i < n; ++i) {
	if (eta(i) == 0.0) return i;
    }
    return n;
}


arma::mat tabulate(const arma::uvec &g, uword j) 
{
    std::map<int, int> table;
    for (uword i = 0; i < g.n_elem; ++i) {
	if (i == j) continue;
	table[g(i)] += 1;
    }

    int n_elem = table.size();
    arma::mat res = arma::mat(2, n_elem, arma::fill::zeros);
    int count = 0;

    for (auto it = table.begin(); it != table.end(); ++it) {
	res(0, count) = it->first;
	res(1, count) = it->second;
	count++;
    }
    return res;
}


uword sample(arma::vec p) 
{
    p /= sum(p);	// normalise
    return find(cumsum(p) > R::runif(0, 1)).eval()(0);
}


// [[Rcpp::export]]
Rcpp::List update_g_marginal(arma::uvec g, arma::vec eta, const arma::vec &y, 
	const arma::vec &x, double sigma, double sigma_0, double a0, 
	Rcpp::Function H0)
{
    g -= 1; // 0-index g so g is compatible with C++ indexing
    uword n = g.n_elem;
    int max_groups = eta.n_elem;
    int cur_group = num_groups(eta);

    for (uword i = 0; i < n; ++i) {
	if (!is_singleton(g, i) && cur_group != max_groups) {
	    double eta_new = Rcpp::as<double>(H0());	// sample from H0
	    double eta_old = eta(g(i));

	    double a = std::min(1.0, a0/(n-1) * 
		    marginal_lik_ratio(eta_new, eta_old, y(i), x(i), sigma, sigma_0));

	    if (a > R::runif(0, 1)) {
		g(i) = cur_group;
		eta(cur_group) = eta_new;
		cur_group += 1;
	    }
	} else {
	    arma::mat table = tabulate(g, i);
	    if (table.n_cols == 1) {
		g(i) = table(0, 0);
	    } else {
		int g_new = table(0, sample(table.row(1).t()));
		int g_old = g(i);
		double a = std::min(1.0, a0/(n-1) * 
		    marginal_lik_ratio(eta(g_new), eta(g_old), y(i), x(i), sigma, sigma_0));
		if(a > R::runif(0, 1)) 
		    g(i) = g_new;
	    }
	}
    }
    
    for (uword i = 0; i < n; ++i) {
	if (!is_singleton(g, i)) {
	    arma::mat table = tabulate(g, i);
	    if (table.n_cols == 1) {
		g(i) = table(0, 0);
	    } else {
		arma::vec probs = table.row(1).t();
		for (uword k = 0; k < probs.n_elem; ++k)
		    probs(k) *= marginal_lik(y(i), x(i), eta(table(0, k)), sigma, sigma_0);

		g(i) = table(0, sample(probs));
	    }
	}
    }
    
    g += 1; // 1-index g, so that g is compatible with R indexing

    return Rcpp::List::create(
	Rcpp::Named("g") = g,
	Rcpp::Named("eta") = eta
    );
}

