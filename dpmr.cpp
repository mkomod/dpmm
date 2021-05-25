#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

using arma::uword;


double dp_plr_lik(double y, double x, double eta, double beta, double sigma)
{
    return R::dnorm(y, eta + beta * x, sigma, 0) ;
}

arma::vec dp_plr_lik(double y, double x, arma::vec eta, double beta, 
	double sigma)
{
    return eta.transform([=] (double e) { 
	return R::dnorm(y, e + beta * x, sigma, 0); 
    });
}


// [[Rcpp::export]]
double dp_plr_log_lik(double y, double x, double eta, double beta, 
	double sigma)
{
    return R::dnorm(y, eta + beta * x, sigma, 1) ;
}


double dp_plr_lik_ratio(double eta_new, double eta_old, double y, double x,
	double beta, double sigma)
{
    double res = dp_plr_lik(y, x, eta_new, beta, sigma) /
	dp_plr_lik(y, x, eta_old, beta, sigma);
    return res;
}


// checks if the jth element in a is a singleton (i.e. unique)
// [[Rcpp::export]]
bool is_singleton(const arma::uvec &a, uword j) 
{
    for (uword i = 0; i < a.n_elem; ++i) {
	if (i == j) continue;
	if (a(i) == a(j)) return false;
    }
    return true;
}

// [[Rcpp::export]]
int num_groups(const arma::vec eta) 
{
    int n = eta.n_elem;
    for (int i = 0; i < n; ++i) {
	if (eta(i) == 0.0) return i;
    }
    return n;
}

arma::mat rdiscrete(int n, arma::vec p, bool as_indicator=false)
{

    if (sum(p) != 1.0)
	Rcpp::stop("probabilities must sum to 1");

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

// [[Rcpp::export]]
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
    p /= sum(p);
    return find(cumsum(p) > R::runif(0, 1)).eval()(0);
}


// [[Rcpp::export]]
Rcpp::List update_g(arma::uvec g, arma::vec eta, const arma::vec &y, 
	const arma::vec &x, double beta, double sigma, double a0, 
	Rcpp::Function H0)
{
    g -= 1; // 0 index g
    uword n = g.n_elem;
    int max_groups = eta.n_elem;
    int cur_group = num_groups(eta);

    for (uword i = 0; i < n; ++i) {
	if (!is_singleton(g, i) && cur_group != max_groups) {
	    double eta_new = Rcpp::as<double>(H0());	// sample from H0
	    double eta_old = eta(g(i));

	    double a = std::min(1.0, a0/(n-1) * 
		    dp_plr_lik_ratio(eta_new, eta_old, y(i), x(i), beta, sigma));

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
		    dp_plr_lik_ratio(eta(g_new), eta(g_old), y(i), x(i), beta, sigma));
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
		    probs(k) *= dp_plr_lik(y(i), x(i), eta(table(0, k)), beta, sigma);

		g(i) = table(0, sample(probs));
	    }
	}
    }
    
    g += 1; // 1 index g, so that g is compatible with R indexing

    return Rcpp::List::create(
	Rcpp::Named("g") = g,
	Rcpp::Named("eta") = eta
    );
}

