#'@examples
#'V=5
#'S=4
#'=3
#'tau_vgb=NULL
#'tau_vgb_0=NULL
#'pi_gs_0=NULL
#'gs="jags"
#'fixed_tau=is.null(tau_vgb)
#'alpha_tau=NULL
#'alpha_epsilon=NULL
#'alpha_bar_epsilon=NULL
#'kappa_rho=NULL
#'alpha_rho=NULL
#'bar_epsilon_1_std=NULL
#'bar_epsilon_1_mean=NULL
#'alpha_pi=.1
#'smc_inits(V,
#'          S,
#'          G,
#'          alpha_tau=.4,
#'          alpha_pi=alpha_pi) 
#'smc_inits(V,
#'          S,
#'          G,
#'          alpha_tau=.4,
#'          alpha_pi=alpha_pi) 
smc_inits <- 
  function(v,
           s,
           g,
           i=1,
           gs="jags",
           tau_vgb_0=NULL,#set initial value
           fixed_tau=FALSE,
           tau_ivgb_0=NULL,#set initial value
           alpha_tau=NULL,
           pi_igs_0=NULL,#set initial value
           pi_gs_0=NULL,#set initial value
           fixed_pi=FALSE,
           alpha_pi=NULL,
           bar_epsilon_1_0=NULL,
           alpha_bar_epsilon=NULL,
           bar_epsilon_1_std=NULL,
           bar_epsilon_1_mean=NULL,
           fixed_bar_epsilon=FALSE,
           epsilon_ba_0=NULL,
           alpha_epsilon=NULL,
           kappa_rho=NULL,
           alpha_rho=NULL) {
    
    
    theta=list()
    sim=if(i>1){sim_tau_pi_epsilon_n_i(v=v,g=g,s=s,i=i,n=0)
      }else{sim_tau_pi_epsilon_n(v=v,g=g,s=s,n=0)}
    
    
    if(is.null(kappa_rho)&is.null(alpha_rho)){
      if(!fixed_bar_epsilon&
         is.null(alpha_bar_epsilon)&
         !is.null(bar_epsilon_1_std)&
         !is.null(bar_epsilon_1_mean)){
        alpha_bar_epsilon=alpha_bar_epsilon_specification(
          bar_epsilon_1_std =bar_epsilon_1_std,
          bar_epsilon_1_mean=bar_epsilon_1_mean)}
      if(!(fixed_bar_epsilon)&!is.null(alpha_bar_epsilon)){
        
        theta=c(theta,        
                if(i>1){list(epsilon_iba=sim$epsilon_iba,
                             bar_epsilon_1_i=3*sim$epsilon_iba[,2,1])},
                if(i==1){list(epsilon_ba=sim$epsilon_ba,
                              bar_epsilon_1=3*sim$epsilon_ba[2,1])})
      }
      if(!is.null(tau_vgb_0)){theta=c(theta,list(tau_vgb=tau_vgb_0))}
      if(!is.null(tau_ivgb_0)){theta=c(theta,list(tau_ivgb=tau_ivgb_0))}
      if(i==1&is.null(tau_vgb_0)&!fixed_tau){theta=c(theta,list(tau_vgb=sim$tau_vgb))}
      if(i>1&is.null(tau_ivgb_0)&!fixed_tau){theta=c(theta,list(tau_ivgb=sim$tau_ivgb))}
      if(!is.null(kappa_rho)){alpha_rho=sample_alpha_rho(kappa_rho);
      theta=c(theta,list(alpha_rho=sample_alpha_rho(kappa_rho)))}}
    if(!is.null(alpha_rho)){
      rep_alpha_rho=rep(alpha_rho,4)
      theta=c(theta,         if(i>1){list(rho_ivga=sim$rho_ivga)},
              if(i==1){list(rho_vgb=sim$rho_vgb)})}
    theta=c(theta,
            if(!is.null(pi_gs_0)){list(pi_gs=pi_gs_0)},
            if(!is.null(pi_igs_0)){list(pi_igs=pi_igs_0)},
            if(i>1&is.null(pi_igs_0)){list(pi_igs=sim$pi_igs)},
            if(i==1&is.null(pi_gs_0)){list(pi_gs=sim$pi_gs)})
    theta
  }