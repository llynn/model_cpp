#include <Rcpp.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
//#include <math.h>
#include <stdlib.h>
//#include <time.h>
#include <R_ext/Utils.h>
using namespace std;

void updatebeta_scalar(vector<double>& xmust, vector<double>& xmuut, vector<int>& xn_s, vector<int>& xn_u, int xI, int xK, int xM,
                      int K1, int t1, double xp_var, double sqrt_var1, double sqrt_var2, int xtt, vector<int>& xgammat, 
                      Rcpp::IntegerMatrix& xd, Rcpp::NumericMatrix& xybar_s, Rcpp::NumericMatrix& xybar_u, Rcpp::NumericMatrix& xys2_s, 
                      Rcpp::NumericMatrix& xys2_u, double xlambda, Rcpp::NumericVector& beta, double xalpha, double xsig_beta1, 
                      double xbeta1, Rcpp::IntegerVector& xAbeta);