#include <Rcpp.h>
#include <R_ext/Utils.h>
#include <stdlib.h>
using namespace std;

void updategamma_indi(vector<int>& xn_s,vector<int>& xn_u, vector<int>& gamma_tt, int xI, int xK, int xM, int xSS, int K1,
                     vector<double>& xalphau, vector<double>& xalphas, double xa, double xb, vector<double>& xmk, double& xIstar, double& xmKstar, 
                     vector<double>& xpp, double xpb1, double xpb2, double xlambda, double xbeta, double xalpha, vector<double>& xmu_s, 
                     vector<double>& xmu_u, Rcpp::IntegerMatrix& xd, Rcpp::NumericMatrix& xybar_s, Rcpp::NumericMatrix& xybar_u, 
                     Rcpp::NumericMatrix& xys2_s, Rcpp::NumericMatrix& xys2_u, Rcpp::IntegerMatrix& xindicator, vector<int>& xAg,
                     double ab, double abI, double bI);