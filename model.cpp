#include <Rcpp.h>
#include <gsl_rng.h>
#include <gsl_randist.h>

#include <stdlib.h>

#include "updatealphau.h"
#include "updategamma_indi.h"
#include "updatealphas.h"
#include "updatemus.h"
#include "updatemuu.h"
#include "updatealpha_scalar.h"
#include "updatebeta_scalar.h"
#include <R_ext/Utils.h>

using namespace std;

RcppExport SEXP model(SEXP T, SEXP I, SEXP K, SEXP M, SEXP ttt, SEXP SS, SEXP alpha_u, SEXP alpha_s, SEXP mu_u, SEXP mu_s, SEXP alpha, 
                      SEXP beta, SEXP gamma, SEXP n_s, SEXP n_u, SEXP varp_u, SEXP lambda_u,SEXP indi, SEXP d, SEXP ybar_s, SEXP ybar_u,
                      SEXP ys2_s, SEXP ys2_u, SEXP a, SEXP b, SEXP lambda, SEXP mk, SEXP Istar, SEXP mKstar, SEXP pp, SEXP pb1,
                      SEXP pb2, SEXP lambda_s,SEXP var_1, SEXP var_2, SEXP p_var, SEXP p_vars, SEXP var_1s, SEXP var_2s, SEXP m_s,SEXP Sigma_s,
                      SEXP p_varu, SEXP var_1u,SEXP var_2u, SEXP m_u, SEXP Sigma_u, SEXP p_vara, SEXP var_1a, SEXP var_2a, SEXP sig_alpha1, SEXP alpha1,
                      SEXP p_varb, SEXP var_1b, SEXP var_2b, SEXP sig_beta1, SEXP beta1, SEXP A_alphau,  SEXP A_alphas, SEXP A_gm, SEXP A_mus, SEXP A_muu, 
                      SEXP A_alpha, SEXP A_beta)
{
  BEGIN_RCPP
  int xT = Rcpp::as<int>(T); 
  int xI = Rcpp::as<int>(I); 
  int xK = Rcpp::as<int>(K); 
  int xtt = Rcpp::as<int>(ttt);
  int xM = Rcpp::as<int>(M); 
  int K1 = xK-1; 
  int xSS = Rcpp::as<int>(SS);
  
  vector<int> xn_s = Rcpp::as<vector<int> >(n_s); 
  vector<int> xn_u = Rcpp::as<vector<int> >(n_u); 
  Rcpp::IntegerMatrix xd(d);
  vector<double> xmk = Rcpp::as<vector<double> >(mk);
  double xIstar = Rcpp::as<double>(Istar);
  double xmKstar = Rcpp::as<double>(mKstar);    
  
  Rcpp::NumericMatrix xybar_s(ybar_s);
  Rcpp::NumericMatrix xybar_u(ybar_u);
  Rcpp::NumericMatrix xys2_s(ys2_s);
  Rcpp::NumericMatrix xys2_u(ys2_u);
  Rcpp::IntegerMatrix xindicator(indi);
  
  // mcmc tuning parameters
  vector<double> sqrt_var = Rcpp::as<vector<double> >(varp_u); //for alpha_u 
  vector<double> xlambda_u =Rcpp::as<vector<double> >(lambda_u);
  double xa = Rcpp::as<double>(a); //for gamma
  double xb = Rcpp::as<double>(b);
  double ab = xa+xb; double abI = ab+xI; double bI = xI+xb; 
  double xlambda = Rcpp::as<double>(lambda); 
  vector<double> xpp = Rcpp::as<vector<double> >(pp); 
  double xpb1 = Rcpp::as<double>(pb1);
  double xpb2 = Rcpp::as<double>(pb2);
  vector<double> xlambda_s = Rcpp::as<vector<double> >(lambda_s);
  vector<double> sqrt_var1 = Rcpp::as<vector<double> >(var_1); //for alpha_s
  vector<double> sqrt_var2 = Rcpp::as<vector<double> >(var_2);
  vector<double> xp_var = Rcpp::as<vector<double> >(p_var); 
  vector<double> sqrt_var1s = Rcpp::as<vector<double> >(var_1s); // for mu_s
  vector<double> sqrt_var2s = Rcpp::as<vector<double> >(var_2s); 
  vector<double> xp_vars = Rcpp::as<vector<double> >(p_vars);
  vector<double> xms = Rcpp::as<vector<double> >(m_s); 
  vector<double> xSigs = Rcpp::as<vector<double> >(Sigma_s);
  vector<double> sqrt_var1u = Rcpp::as<vector<double> >(var_1u); // for mu_u
  vector<double> sqrt_var2u = Rcpp::as<vector<double> >(var_2u);
  vector<double> xp_varu = Rcpp::as<vector<double> >(p_varu);
  vector<double> xmu = Rcpp::as<vector<double> >(m_u); 
  vector<double> xSigu = Rcpp::as<vector<double> >(Sigma_u);
  double xalpha1 = Rcpp::as<double>(alpha1); //for alpha
  double xsig_alpha1 = Rcpp::as<double>(sig_alpha1);
  double sqrt_var1a = Rcpp::as<double>(var_1a);
  double sqrt_var2a = Rcpp::as<double>(var_2a); double xp_vara = Rcpp::as<double>(p_vara);
  double sqrt_var1b = Rcpp::as<double>(var_1b); // for beta
  double sqrt_var2b = Rcpp::as<double>(var_2b);
  double xp_varb = Rcpp::as<double>(p_varb);
  double xbeta1 = Rcpp::as<double>(beta1);
  double xsig_beta1 = Rcpp::as<double>(sig_beta1);


  // Initialize model parameters to store mcmc traces
  Rcpp::NumericMatrix alphaut(alpha_u); 
  Rcpp::IntegerMatrix Aalphau(A_alphau);
  Rcpp::NumericMatrix alphast(alpha_s);
  Rcpp::IntegerMatrix Aalphas(A_alphas);
  Rcpp::NumericMatrix muut(mu_u);  
  Rcpp::IntegerMatrix Amuu(A_muu);
  Rcpp::NumericMatrix must(mu_s); 
  Rcpp::IntegerMatrix Amus(A_mus);
  Rcpp::NumericVector alphat(alpha); 
  Rcpp::IntegerVector Aalpha(A_alpha);
  Rcpp::NumericVector betat(beta); 
  Rcpp::IntegerVector Abeta(A_beta);
  Rcpp::IntegerMatrix gammat(gamma);
  Rcpp::IntegerMatrix Agamma(A_gm);
  
  vector<double> alphau_t1(xK);
  vector<double> alphas_t1(xK);
  vector<int> gamma_t1(xI*xK);
  vector<int> Agm_t1(xI);
  vector<int> Aalp_t1(xK);
  vector<int> Aalps_t1(xK);
  vector<int> Amus_t1(xM); vector<int> Amuu_t1(xM);
  vector<double> mus_t1(xM); vector<double> muu_t1(xM); 
  double log1 = 0.0;  double log2 = 0.0; double delF = 0.0;  
  Rcpp::RNGScope scope;
  int t1 = 0;
  
  for (int tt=1; tt<xT; tt++){
    t1 = tt-1;
    for ( int k=0; k<xK; k++){
      alphau_t1[k] = alphaut(t1,k); 
      Aalp_t1[k] = 0; 
      alphas_t1[k] = alphast(t1,k);
    }
    for (int m=0; m<xM; m++){
      mus_t1[m] = must(t1,m);
      muu_t1[m] = muut(t1,m);
    }
    for ( int ik=0; ik<(xI*xK); ik++) {gamma_t1[ik] = gammat(ik,t1);}
    // update alpha_u
    updatealphau(alphau_t1, xn_s, xn_u, xI, xK, xlambda_u, sqrt_var, xtt, gamma_t1, Aalp_t1);
    for (int k=0; k<xK; k++){
      alphaut(tt,k) = alphau_t1[k]; //alphau_t1 is updated
      Aalphau(k,tt) = Aalp_t1[k];
    }
    
    // update gamma
    updategamma_indi(xn_s,xn_u, gamma_t1, xI,  xK,  xM,  xSS, K1, alphau_t1, alphas_t1, xa, xb,  xmk,  xIstar, xmKstar, 
                     xpp,  xpb1,  xpb2, xlambda, betat[t1], alphat[t1], mus_t1, muu_t1, xd,xybar_s, xybar_u, 
                     xys2_s,  xys2_u, xindicator,  Agm_t1,ab, abI, bI);
    for (int ik=0; ik<(xI*xK); ik++) {
      gammat(ik,tt) = gamma_t1[ik]; //gamma_t1 is updated
    }
    
    // update alpha_s
    updatealphas(alphas_t1, xn_s, xK, xI,  xlambda_s, gamma_t1, sqrt_var1, sqrt_var2, xp_var, xtt, Aalps_t1);
    for (int k=0; k<xK; k++){
      alphast(tt,k) = alphas_t1[k]; //alphas_t1 is updated
      Aalphas(k,tt) = Aalps_t1[k];
    }
    
    //update mu_s
    updatemus(mus_t1, muu_t1, xn_s, xn_u, xI, xK,  xM, K1, xp_vars, sqrt_var1s, sqrt_var2s, xtt, gamma_t1, xd, xybar_s, 
              xybar_u, xys2_s, xys2_u, xlambda, betat[t1], alphat[t1], xms,xSigs , Amus_t1);
   for (int m=0; m<xM; m++){
     must(tt,m) = mus_t1[m];
     Amus(m,tt) = Amus_t1[m];
   }
   
   // update mu_u
   updatemuu(mus_t1,muu_t1, xn_s, xn_u, xI, xK, xM, K1, xp_varu, sqrt_var1u, sqrt_var2u, xtt, gamma_t1, xd, 
             xybar_s, xybar_u, xys2_s, xys2_u, xlambda, betat[t1], alphat[t1], xmu, xSigu, Amuu_t1);
   for (int m=0; m<xM; m++){
     muut(tt,m) = muu_t1[m];
     Amuu(m,tt) = Amuu_t1[m];
   }
    
   // update alpha 
   updatealpha_scalar(mus_t1, muu_t1, xn_s, xn_u, xI, xK, xM,  K1,  t1, xp_vara, sqrt_var1a, sqrt_var2a, xtt, gamma_t1,
                      xd, xybar_s, xybar_u, xys2_s, xys2_u, xlambda, betat[t1], alphat, xsig_alpha1, xalpha1, Aalpha);
                      

  
  // update beta
  updatebeta_scalar(mus_t1, muu_t1, xn_s, xn_u, xI, xK, xM, K1, t1, xp_varb, sqrt_var1b, sqrt_var2b, xtt, gamma_t1, 
                    xd, xybar_s, xybar_u, xys2_s, xys2_u, xlambda, betat, alphat[(t1+1)], xsig_beta1, xbeta1, Abeta);
  }
  
END_RCPP
  
}
  
  
  
  