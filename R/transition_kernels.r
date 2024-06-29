
list(
temper_bar_epsilon=FALSE,
temper_n_vsa=FALSE,
temper_alpha_tau=FALSE)
#'@examples
#'
#'fixed_bar_epsilon=FALSE
#'constrained_epsilon_matrix=TRUE
#'block_tau=FALSE
#'fixed_tau=FALSE
#'relax_tau=TRUE
#'relax_rho=FALSE
#'n=1000
#'v=50
#'s=3
#'g=6
#'rho_ivga
#'theta=sim_tau_pi_epsilon_n_i(i=4,v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)
#'obs=sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)
#'m_ivsabg=sampler_m_ivsabg(theta$chi_ivsabg,n_vsa)
#'(kernel_f(
#'    fixed_bar_epsilon=FALSE,constrained_epsilon_matrix=TRUE,block_tau=FALSE,fixed_tau=FALSE,relax_tau=TRUE,relax_rho=FALSE))(
#'          theta=theta,
#'          fixed=list(n_vsa=n_vsa,rep_alpha_pi=rep(.1,g),rep_alpha_tau=rep(.1,4),alpha_bar_epsilon=c(1,100)))
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
    fixed){
            allvar=c(theta,fixed);
            thetanew=list();
                             ",
    if(relax_rho){
      "
      chi_ivsag=einsum::einsum('ivga,igs->ivsag',allvar$rho_ivga,allvar$pi_igs);
      m_ivsag=sampler_m_ivsag(n_vsa=allvar$n_vsa,chi_ivsag=chi_ivsag);
      thetanew$rho_ivga<-
        rho_ivga<-
        sampler_tilde_rho_ivga(m_ivga=einsum::einsum('ivsag->ivga',m_ivsag),
                               rep_alpha_rho=allvar$rep_alpha_rho);
      thetanew$pi_igs<-
        pi_igs<-
        sampler_pi_igs(m_igs=einsum::einsum('ivsag->igs',m_ivsag),
                       rep_alpha_pi=allvar$rep_alpha_pi);"
    },
    if(!relax_rho){
      "
      chi_ivsabg=einsum::einsum('ivgb,iba,igs->ivsabg',allvar$tau_ivgb,allvar$epsilon_iba,allvar$pi_igs);
      m_ivsabg=sampler_m_ivsabg(n_vsa=n_vsa,chi_ivsabg=chi_ivsabg);
       thetanew$pi_igs<-
          pi_igs<-
          sampler_pi_igs(m_igs=einsum::einsum('ivsabg->igs',m_ivsabg),
                         rep_alpha_pi=allvar$rep_alpha_pi);"},
    if(!relax_rho&relax_tau&!fixed_tau){
      "thetanew$tau_ivgb<-
        tau_ivgb<-
        sampler_tilde_tau_ivgb(m_ivgb=einsum::einsum('ivsabg->ivgb',m_ivsabg),
                               rep_alpha_tau=allvar$rep_alpha_tau);"
      },
    if(!relax_rho&!relax_tau&!fixed_tau&!block_tau){
      "tau_ivgb<-
          thetanew$tau_ivgb<-
          sampler_tau_ivgb(theta$tau_ivgb,
                        pi_igs,
                        allvar$epsilon_iba,
                        n_vsa,
                        block_tau=FALSE,
                        v=dim(theta$tau_ivgb)[2],
                        #s=dim(theta$pi_gs)[2],
                        g=dim(theta$pi_igs)[2],
                        g_neq_g=g_neq_g_f(g),
                        alpha_tau=0,
                        m_ivgb=einsum::einsum('ivsabg->ivgb',m_ivsabg),
                           rep_alpha_tau=allvar$rep_alpha_tau);"
    },
    if(!relax_rho&!relax_tau&!fixed_tau&block_tau){
      "tau_ivgb<-
          thetanew$tau_ivgb<-
          sampler_tau_ivgb(theta$tau_ivgb,
                        pi_igs,
                        epsilon_iba,
                        n_vsa,
                        block_tau=TRUE,
                        v=dim(theta$tau_ivgb)[2],
                        #s=dim(theta$pi_gs)[2],
                        g=dim(theta$pi_igs)[2],
                        g_neq_g=g_neq_g_f(g),
                        alpha_tau=0,
                        m_ivgb=einsum::einsum('ivsabg->ivgb',m_ivsabg),
                           rep_alpha_tau=allvar$rep_alpha_tau);"
    },
    if(!relax_rho&!fixed_bar_epsilon&constrained_epsilon_matrix){
      "m_ianeqb=einsum::einsum('ivsabg,ab->i',m_ivsabg,a_neq_b);
       m_iaa=einsum::einsum('ivsaag->i',m_ivsabg);
       bar_epsilon_1i=sampler_bar_epsilon_i1(m_ianeqb=m_ianeqb,m_iaa=m_iaa,alpha_bar_epsilon=allvar$alpha_bar_epsilon);
      allvar$epsilon_iba<-thetanew$epsilon_iba<-epsilon_iba<-plyr::aaply(bar_epsilon_1i,1,epsilon_ba_f);
      "
    },
    if(!relax_rho&!fixed_bar_epsilon&!constrained_epsilon_matrix){
      "m_iba=einsum::einsum('ivsabg->iba',m_ivsabg);
      thetanew$epsilon_iba<-epsilon_iba<-sampler_epsilon_iba(m_iba=m_iba,rep_alpha_epsilon=allvar$rep_alpha_epsilon);
      "
    },
    if(!relax_rho&!fixed_bar_epsilon&!constrained_epsilon_matrix){
      "m_iba=einsum::einsum('ivsabg->iba',m_ivsabg);
      allvar$epsilon_iba<-epsilon_iba<-sampler_epsilon_iba(m_iba=m_iba,rep_alpha_epsilon=allvar$rep_alpha_epsilon);
      "
    },
    "thetanew}")))
    
    
  }
#'@description SMC kernel
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha_pi=alpha_pi)|>attach()

desman_kernel<-function(n_vsa,tau_vgb,pi_gs,epsilon_ba,rep_alpha_pi,delta){
  xi_vsabg=sampler_m(tau_vgb=tau_vgb,
                      pi_gs=pi_gs,
                      epsilon_ba = epsilon_ba)
  nu_vsab<-m_vsab_from_m(xi=xi_vsabg)
  mu_vsab<-m_vsag_from_m(xi=xi_vsabg)
  pi_gs<-sampler_pi(mu_vsag = mu_vsab,alpha_pi = rep_alpha_pi)
  w1<-smc_b_prime(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
  tau_vgb<-sampler_tau(n_vsa = n_vsa,tau_vgb = tau_vgb,pi_gs = pi_gs,epsilon_ba = epsilon_ba)
  w2<-smc_b_prime(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
  epsilon_ba<-sampler_epsilon_star(nu_vsab,delta)
  w3<-smc_b_prime(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
  list(tau_vgb=tau_vgb,pi_gs=pi_gs,epsilon_ba=epsilon_ba,w1=w1,w2=w2,w3=w3)
}
