#include <Rcpp.h>
#include <gsl_rng.h>
#include <gsl_randist.h>
#include <stdlib.h>


#include <R_ext/Utils.h>
using namespace std;

void updatealphas(vector<double>& xalphast,vector<int>& xn_s, int xK, int xI, vector<double>& xlambda_s, vector<int>& xgammat, 
                  vector<double>& sqrt_var1,vector<double>& sqrt_var2, vector<double>& xp_var, int xtt, vector<int>& xAalphas);