#include <RcppArmadillo.h>
#include <Rmath.h>
#include "utils.h" 
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace Rcpp;

List ddloglik_exact(const arma::mat& Z,
                    const arma::vec& delta,
                    const arma::vec& time,
                    const arma::vec& beta,
                    const arma::vec& n_each_stratum,
                    const double comb_max) {

    const int p = Z.n_cols;
    const int S = n_each_stratum.n_elem;
    
    arma::uvec ind_start(S);
    ind_start(0) = 0u;
    for (int s = 1; s < S; ++s) {
        ind_start(s) = ind_start(s - 1) + (unsigned int)n_each_stratum(s - 1);
    }
    
    arma::vec loglik(1, arma::fill::zeros);
    arma::vec L1(p, arma::fill::zeros);
    arma::mat L2(p, p, arma::fill::zeros);
    
    for (int s = 0; s < S; ++s) {
        int start = (int)ind_start(s);
        int len = (int)n_each_stratum(s);
        int end = start + len - 1;
        
        arma::mat Z_s = Z.rows(start, end);
        arma::vec delta_s = delta.subvec(start, end);
        arma::vec time_s = time.subvec(start, end);
        
        int i = 0;
        while (i < len) {
            double tk = time_s(i);
            int b = i;
            while (i < len && time_s(i) == tk) ++i;
            int e = i - 1;
            
            int d_k = 0;
            for (int r = b; r <= e; ++r) if (delta_s(r) == 1.0) ++d_k;
            if (d_k == 0) continue;
            
            arma::rowvec w_k(p, arma::fill::zeros);
            for (int r = b; r <= e; ++r) if (delta_s(r) == 1.0) w_k += Z_s.row(r);
            
            int n_risk = len - b;
            arma::mat Z_risk = Z_s.rows(b, len - 1);
            
            double ncomb_dbl = R::choose((double)n_risk, (double)d_k);
            if (ncomb_dbl > comb_max) {
                stop("Exact enumeration too large. Increase comb_max or use Efron/Breslow.");
            }
            
            arma::umat combs = combn_index(n_risk, d_k);
            const int n_comb = combs.n_cols;
            
            double eta_max = -std::numeric_limits<double>::infinity();
            for (int j = 0; j < n_comb; ++j) {
                arma::uvec idx = combs.col(j);
                arma::rowvec w = arma::sum(Z_risk.rows(idx), 0);
                double eta = arma::dot(w, beta);
                if (eta > eta_max) eta_max = eta;
            }
            
            double S0 = 0.0;
            arma::rowvec S1(p, arma::fill::zeros);
            arma::mat S2(p, p, arma::fill::zeros);
            
            for (int j = 0; j < n_comb; ++j) {
                arma::uvec idx = combs.col(j);
                arma::rowvec w = arma::sum(Z_risk.rows(idx), 0);
                double eta = arma::dot(w, beta);
                double wgt = std::exp(eta - eta_max);
                
                S0 += wgt;
                S1 += wgt * w;
                S2 += wgt * (w.t() * w);
            }
            
            loglik(0) += arma::dot(w_k, beta) - (eta_max + std::log(S0));
            arma::rowvec score_k = w_k - (S1 / S0);
            L1 += score_k.t();
            arma::mat mean_outer = (S1.t() * S1) / (S0 * S0);
            arma::mat hess_k = (S2 / S0) - mean_outer;
            L2 += hess_k;
        }
    }
    
    return List::create(Named("loglik") = loglik, Named("L1") = L1, Named("L2") = L2);
}


List ddloglik_approx(const arma::mat& Z,
                     const arma::vec& delta,
                     const arma::vec& time,
                     const arma::vec& beta,
                     const arma::vec& n_each_stratum,
                     int method) {
    
    const int p = Z.n_cols;
    const int S = n_each_stratum.n_elem;
    
    arma::uvec ind_start(S);
    ind_start(0) = 0u;
    for (int s = 1; s < S; ++s) {
        ind_start(s) = ind_start(s - 1) + (unsigned int)n_each_stratum(s - 1);
    }
    
    arma::vec loglik(1, arma::fill::zeros);
    arma::vec L1(p, arma::fill::zeros);
    arma::mat L2(p, p, arma::fill::zeros);
    
    for (int s = 0; s < S; ++s) {
        int start = (int)ind_start(s);
        int len = (int)n_each_stratum(s);
        int end = start + len - 1;
        
        arma::mat Z_s = Z.rows(start, end);
        arma::vec delta_s = delta.subvec(start, end);
        arma::vec time_s = time.subvec(start, end);
        
        arma::vec theta_s = Z_s * beta;
        arma::vec exp_theta_s = arma::exp(theta_s);
        arma::vec S0_s = rev_cumsum(exp_theta_s);
        
        arma::mat S1_s(len, p, arma::fill::zeros);
        for (int j = 0; j < p; ++j) {
            arma::vec pre = Z_s.col(j) % exp_theta_s;
            S1_s.col(j) = rev_cumsum(pre);
        }
        
        std::vector<int> b_idx;
        std::vector<int> d_vals;
        std::vector<arma::rowvec> wk_list;
        std::vector<double> E0_list;
        std::vector<arma::rowvec> E1_list;
        std::vector<arma::mat> E2_list;
        
        int i = 0;
        while (i < len) {
            double tk = time_s(i);
            int b = i;
            while (i < len && time_s(i) == tk) ++i;
            int e = i - 1;
            
            int d_k = 0;
            arma::rowvec w_k(p, arma::fill::zeros);
            double E0 = 0.0;
            arma::rowvec E1(p, arma::fill::zeros);
            arma::mat E2(p, p, arma::fill::zeros);
            
            for (int r = b; r <= e; ++r) {
                if (delta_s(r) == 1.0) {
                    ++d_k;
                    w_k += Z_s.row(r);
                    if (method == 1) { 
                        double ef = exp_theta_s(r);
                        E0 += ef;
                        E1 += ef * Z_s.row(r);
                        arma::rowvec zr = Z_s.row(r);
                        E2 += ef * (zr.t() * zr);
                    }
                }
            }
            if (d_k == 0) continue;
            
            b_idx.push_back(b);
            d_vals.push_back(d_k);
            wk_list.push_back(w_k);
            if (method == 1) {
                E0_list.push_back(E0);
                E1_list.push_back(E1);
                E2_list.push_back(E2);
            }
        }
        
        const int K = (int)b_idx.size();
        for (int k = 0; k < K; ++k) {
            int b = b_idx[k];
            int d = d_vals[k];
            const arma::rowvec& w_k = wk_list[k];
            double S0_base = S0_s(b);
            arma::rowvec S1_base = S1_s.row(b);
            
            if (method == 0) { 
                loglik(0) += arma::dot(w_k, beta) - d * std::log(S0_base);
                arma::rowvec score_k = w_k - d * (S1_base / S0_base);
                L1 += score_k.t();
            } else { 
                double E0 = E0_list[k];
                arma::rowvec E1 = E1_list[k];
                double log_denom = 0.0;
                arma::rowvec adj(p, arma::fill::zeros);
                for (int r = 0; r < d; ++r) {
                    double a = (double)r / (double)d;
                    double S0r = S0_base - a * E0;
                    log_denom += std::log(S0r);
                    arma::rowvec S1r = S1_base - a * E1;
                    adj += S1r / S0r;
                }
                loglik(0) += arma::dot(w_k, beta) - log_denom;
                L1 += (w_k - adj).t();
            }
        }
        
        for (int a = 0; a < p; ++a) {
            for (int bcol = 0; bcol < p; ++bcol) {
                arma::vec S2_vec = rev_cumsum(Z_s.col(a) % exp_theta_s % Z_s.col(bcol));
                double acc = 0.0;
                for (int k = 0; k < K; ++k) {
                    int b = b_idx[k];
                    int d = d_vals[k];
                    double S0_base = S0_s(b);
                    double S2_base = S2_vec(b);
                    double S1a_base = S1_s(b, a);
                    double S1b_base = S1_s(b, bcol);
                    
                    if (method == 0) {
                        double Vij = (S2_base / S0_base) - (S1a_base * S1b_base) / (S0_base * S0_base);
                        acc += d * Vij;
                    } else {
                        double E0 = E0_list[k];
                        double E2ab = E2_list[k](a, bcol);
                        double E1a = E1_list[k](a);
                        double E1b = E1_list[k](bcol);
                        for (int r = 0; r < d; ++r) {
                            double tfrac = (double)r / (double)d;
                            double S0r = S0_base - tfrac * E0;
                            double S2r = S2_base - tfrac * E2ab;
                            double S1ar = S1a_base - tfrac * E1a;
                            double S1br = S1b_base - tfrac * E1b;
                            acc += (S2r / S0r) - (S1ar * S1br) / (S0r * S0r);
                        }
                    }
                }
                L2(a, bcol) += acc;
            }
        }
    }
    
    return List::create(Named("loglik") = loglik, Named("L1") = L1, Named("L2") = L2);
}


// [[Rcpp::export]]
List Cox_NR(const arma::mat& Z,
            const arma::vec& delta,
            const arma::vec& time,
            const arma::vec& n_each_stratum,
            arma::vec beta,
            int ties_method,
            double tol,
            int max_iter,
            double comb_max = 1e7) {
    
    arma::vec update(beta.n_elem);
    arma::vec L1;
    arma::mat L2;
    List diff;
    
    for (int iter = 0; iter < max_iter; ++iter) {
        if (ties_method == 2) {
            diff = ddloglik_exact(Z, delta, time, beta, n_each_stratum, comb_max);
        } else {
            diff = ddloglik_approx(Z, delta, time, beta, n_each_stratum, ties_method);
        }
        
        L1 = as<arma::vec>(diff["L1"]);
        L2 = as<arma::mat>(diff["L2"]);
        
        update = arma::solve(L2, L1);
        beta += update;
        
        if (arma::max(arma::abs(update)) < tol) break;
    }
    
    if (ties_method == 2) {
        diff = ddloglik_exact(Z, delta, time, beta, n_each_stratum, comb_max);
    } else {
        diff = ddloglik_approx(Z, delta, time, beta, n_each_stratum, ties_method);
    }
    
    double loglik = as<arma::vec>(diff["loglik"])(0);
    
    return List::create(Named("beta") = beta, Named("loglik") = loglik);
}