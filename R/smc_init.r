sample_inits_generic_f <- 
  function(V,
           S,
           G,
           tau_vgb_0=NULL,
           pi_gs_0=NULL,
           gs="jags",
           fixed_tau_vgb=FALSE,
           alpha_tau=NULL,
           alpha_epsilon=NULL,
           alpha_bar_epsilon=NULL,
           kappa_rho=NULL,
           bar_epsilon_1_std=NULL,
           bar_epsilon_1_mean=NULL,
           alpha_pi=NULL) {
    
    v = dim(n_vsa)[1]
    s = dim(n_vsa)[2]
    
    inits=list()
    if(is.null(kappa_rho)){
    if(is.null(alpha_bar_epsilon)&
       !is.null(bar_epsilon_1_std)&
       !is.null(bar_epsilon_1_mean)){
      alpha_bar_epsilon=alpha_bar_epsilon_specification(
        bar_epsilon_1_std =bar_epsilon_1_std,
                                                  bar_epsilon_1_mean=bar_epsilon_1_mean)}
    if(!(fixed_bar_epsilon_1)&!is.null(alpha_bar_epsilon)){
      inits=c(inits,list(bar_epsilon_1=sampler_bar_epsilon_1_0(alpha_bar_epsilon)))
    }
    inits=c(inits,list(tau_vgb=if(!is.null(tau_vgb)){tau_vgb}else{sim_tau_vgb(v,g)}))}
    if(!is.null(kappa_rho)){inits=c(inits,list(alpha_rho=sample_alpha_rho(kappa_rho),
                                               rho=sample_rho(v,G,alpha_rho)))}
    
    inits=c(inits,list(pi_gs=sim_pi_gs(g,s,alpha_pi)))
    inits}