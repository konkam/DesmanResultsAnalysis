
list(
temper_bar_epsilon=FALSE,
temper_n_vsa=FALSE,
temper_alpha_tau=FALSE)

kernel_f<-function(
        fixed_bar_epsilon=fixed_bar_epsilon,
        constrained_epsilon_matrix=constrained_epsilon_matrix,
        block_tau=block_tau,
        fixed_tau=fixed_tau,
        relax_tau=relax_tau,
        relax_rho=relax_rho){
      eval(parse(text=paste0("
  function(
    theta,
    lambda,
    obs,
    param){",
    if(relax_rho){
      "chi_ivsag=1
      m_ivsag=sampler_m_ivsag(obs$n_vsa,chi_ivsag)
      theta$rho_ivga=sampler_rho_ivga(m_ivsag,alpha_rho)
      m_igs=plyr::aaply(m_ivsag,c(1,5,2),sum)
      theta_pi_igs=sampler_theta_pi_igs(,param$alpha_pi)
      theta$"
    },
    if(!relax_rho){
      "chi_ivsabg=1
      m=sampler_m_ivsabg(obs$n_vsa,chi_ivsabg)
      theta$tau_ivga=sampler_tau_ivga()
      
      theta$"
    },
    if(!relax_tau){
      "chi_ivsabg=1
      m=sampler_m_ivsabg(obs$n_vsa,chi_ivsabg)
      theta$tau_ivga=sampler_tau_ivga()
      
      theta$"
    },
    "theta}")))
    
    
  }
#'@description SMC kernel
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha_pi=alpha_pi)|>attach()




desman_kernel<-function(n_vsa,tau_vgb,pi_gs,epsilon_ba,rep_alpha_pi,delta){
  xi_vsabg=sampler_xi(tau_vgb=tau_vgb,
                      pi_gs=pi_gs,
                      epsilon_ba = epsilon_ba)
  nu_vsab<-nu_from_xi(xi=xi_vsabg)
  mu_vsab<-mu_from_xi(xi=xi_vsabg)
  pi_gs<-sampler_pi(mu_vsag = mu_vsab,alpha_pi = rep_alpha_pi)
  w1<-smc_b_prime(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
  tau_vgb<-sampler_tau(n_vsa = n_vsa,tau_vgb = tau_vgb,pi_gs = pi_gs,epsilon_ba = epsilon_ba)
  w2<-smc_b_prime(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
  epsilon_ba<-sampler_epsilon_star(nu_vsab,delta)
  w3<-smc_b_prime(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
  list(tau_vgb=tau_vgb,pi_gs=pi_gs,epsilon_ba=epsilon_ba,w1=w1,w2=w2,w3=w3)
}

smc_kernel_i<-function(i,g,n_vsa,tau_ivga,pi_igs,epsilon_iba,rep_alpha_pi,delta,g_neq_g =g_neq_g_f(g)){
  xi_ivsabg=sampler_xi_i(i =i, tau_ivga=tau_ivga,
                         pi_igs=pi_igs,
                         epsilon_iba = epsilon_iba)
  nu_ivsab=nu_i_from_xi_i(xi=xi_ivsabg)
  mu_ivsab=mu_i_from_xi_i(xi=xi_ivsabg)
  pi_igs=sampler_pi_i(mu_ivsag = mu_ivsab,rep_alpha_pi = rep_alpha_pi)
  tau_ivga=sampler_tau_i(n_vsa = n_vsa,tau_ivga = tau_ivga,pi_igs = pi_igs,epsilon_iba = epsilon_iba,g_neq_g=g_neq_g )
  list(tau_ivga=tau_ivga,pi_igs=pi_igs,epsilon_iba=epsilon_iba)
}


