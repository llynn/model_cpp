####################blackrhino CD4 #############
set.seed(680)
source("repmat.R")
source("realdata_rv144.R")
source("init_realdata_rv144.R")


dyn.load("model.so")
T = 50;
ttt=2000;
n_ss = as.integer(n_s); n_uu = as.integer(n_u) 
ybars = apply(ybar_s, 3, function(x) as.vector(x))
ybaru = apply(ybar_u, 3, function(x) as.vector(x))
ys2s = apply(ys2_s, 3, function(x) as.vector(x))
ys2u = apply(ys2_u, 3, function(x) as.vector(x))

gammat = apply(gamma,3, function(x) as.integer(x))

ptm <- proc.time()
result<-.Call("model",T=T,I=I, K=K, M=M, ttt=ttt, SS=1,alpha_u=alpha_u, alpha_s=alpha_s, mu_u=mu_u, mu_s=mu_s, alpha=alpha,
                      beta=beta, gamma=gammat,n_s=n_ss, n_u=n_uu, varp_u=varp_u, lambda_u=lambda_u, indi=indi, d=d, ybar_s = ybars, ybar_u=ybaru, 
                      ys2_s = ys2s, ys2_u = ys2u, a=a, b=b,lambda=lambda,mk = mk, Istar = Istar, mKstar = mKstar, pp = pp, pb1 = pb1, 
                      pb2 = pb2,lambda_s = lambda_s, var_1 = varp_s1,var_2 = varp_s2,p_var = pvar_s, p_vars = pvar_mus,var_1s = sqrt_mus1, var_2s=sqrt_mus2,m_s=m_s,Sigma_s =Sigma_s,
                      p_varu=pvar_muu,var_1u=sqrt_muu1,var_2u=sqrt_muu2,m_u=m_u, Sigma_u=Sigma_u,p_vara=pvar_alpha, var_1a=sqrt_alpha1,var_2a=sqrt_alpha2,sig_alpha1=sig_alpha1,alpha1=alpha1,
                      p_varb=pvar_beta,var_1b=sqrt_beta1,var_2b=sqrt_beta2,sig_beta1 = sig_beta1, beta1= beta1, A_alphau=A_alphau,A_alphas=A_alphas, A_gm=A_gm, A_mus=A_mus, A_muu=A_muu,  
                      A_alpha=A_alpha,A_beta=A_beta)
proc.time() - ptm
