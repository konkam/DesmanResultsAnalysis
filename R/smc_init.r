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
           gs="jags",
           tau_vgb_0=NULL,#set initial value
           fixed_tau=FALSE,
           alpha_tau=NULL,
           pi_gs_0=NULL,#set initial value
           fixed_pi=FALSE,
           alpha_pi=NULL,
           bar_epsilon_1_0=NULL,
           alpha_bar_epsilon=NULL,
           bar_epsilon_1_std=NULL,
           bar_epsilon_1_mean=NULL,
           fixed_bar_epsilon_1=FALSE,
           epsilon_ba_0=NULL,
           alpha_epsilon=NULL,
           kappa_rho=NULL,
           alpha_rho=NULL) {
    
    inits=list()
    if(is.null(kappa_rho)&is.null(alpha_rho)){
    if(!fixed_bar_epsilon_1&
       is.null(alpha_bar_epsilon)&
       !is.null(bar_epsilon_1_std)&
       !is.null(bar_epsilon_1_mean)){
      alpha_bar_epsilon=alpha_bar_epsilon_specification(
        bar_epsilon_1_std =bar_epsilon_1_std,bar_epsilon_1_mean=bar_epsilon_1_mean)}
    if(!(fixed_bar_epsilon_1)&!is.null(alpha_bar_epsilon)){
      inits=c(inits,list(bar_epsilon_1=sampler_bar_epsilon_1_0(alpha_bar_epsilon)))
    }
    inits=c(inits,list(tau_vgb=if(!is.null(tau_vgb)){tau_vgb}else{sim_tau_vgb(v,g)}))}
    if(!is.null(kappa_rho)){inits=c(inits,list(alpha_rho=sample_alpha_rho(kappa_rho),
                                               rho=sample_rho(v,g,alpha_rho)))}
    
    if(!is.null(alpha_rho)){inits=c(inits,list(rho=sample_rho(v,g,alpha_rho)))}
    inits=c(inits,list(pi_gs=sim_pi_gs(g,s,alpha_pi)))
    inits}
