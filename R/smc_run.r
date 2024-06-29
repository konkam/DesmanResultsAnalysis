#'@examples
#'gs="custom"
#'tau_pi_n <- sim_tau_pi_epsilon_n(v = 50, g = 5, s = 3, n = 1000, alpha_pi = 1)
#'n_vsa = tau_pi_n$n_vsa
#'tau_vgb = tau_pi_n$tau_vgb
#'G=if(!is.null(tau_vgb)){dim(tau_vgb)[2]}else{5}
#'bar_epsilon_1=NULL
#'tau_vgb=NULL
#'alpha_tau=NULL
#'alpha_epsilon=NULL
#'alpha_bar_epsilon=c(1,10)
#'kappa_rho=NULL
#'kappa_pi=c(1,10)
#'bar_epsilon_1_std=NULL
#'bar_epsilon_1_mean=NULL
#'alpha_pi=1
#'fixed_bar_epsilon=!is.null(bar_epsilon_1)
#'constrained_epsilon_matrix=TRUE
#'block_tau=TRUE
#'relax_tau=!is.null(alpha_tau)
#'relax_rho=!is.null(kappa_rho)
#'n_chains = 2
#'for(gs in c("jags","custom")){
#'assign(paste0("X",gs),
#'smc_run(n_vsa,
#'tau_vgb=tau_vgb,
#'G=G,
#'gs=gs,
#'adapt=0,
#'sample=10)
#')}
#'X|>mcmc_output_df("bar_epsilon")|>View()
    
smc_run<-
  function(n_vsa,
           gs="jags",
           tau_ivgb_0=NULL,
           tau_vgb_0=NULL,
           pi_igs_0=NULL,#set initial value
           pi_gs_0=NULL,#set initial value
           tau_vgb=NULL,
           G=if(!is.null(tau_vgb)){dim(tau_vgb)[2]}else{5},
           block_tau=TRUE,
           alpha_tau=NULL,
           alpha_epsilon=NULL,
           bar_epsilon_1_std=NULL,
           bar_epsilon_1_mean=NULL,
           alpha_bar_epsilon=c(1,10),
           bar_epsilon_1=NULL,
           kappa_rho=NULL,
           alpha_rho=NULL,
           alpha_pi=.1,
           n_chains = 2,
           mcmc=FALSE,
           inits=NULL,
           t_min=1,
           t_max=30,
           ess_min=NULL,
           smc_kernel=desman_kernel,
           trace_all=TRUE,
           .update_lambda=update_lambda,
           n_vsa_df=NULL,
           tempering_n=FALSE,
           bar_epsilon_1_tempering=FALSE,
           alpha_tau_tempering=FALSE,
           alpha_tau_seq=NULL,
           bar_epsilon_1_seq=NULL,
           ...) {
    
    fixed_bar_epsilon=!is.null(bar_epsilon_1)
    fixed_tau=!is.null(tau_vgb)
    relax_tau=!is.null(alpha_tau)
    relax_rho=!is.null(kappa_rho)|!is.null(alpha_rho)
    fixed_alpha_rho=!is.null(alpha_rho)
    constrained_epsilon_matrix=is.null(alpha_epsilon)
    
    model_string <-
      model_string_f(
        gs=gs,
        fixed_bar_epsilon=fixed_bar_epsilon,
        constrained_epsilon_matrix=constrained_epsilon_matrix,
        block_tau=block_tau,
        fixed_tau=fixed_tau,
        relax_tau=relax_tau,
        relax_rho=relax_rho)
    
    observations_and_constants <- 
      smc_fixed_f(n_vsa=n_vsa,
                  gs=gs,
                  G=G,
                  bar_epsilon_1=bar_epsilon_1,
                  tau_vgb=tau_vgb,
                  alpha_tau=alpha_tau,
                  alpha_epsilon=alpha_epsilon,
                  alpha_bar_epsilon=alpha_bar_epsilon,
                  kappa_rho=kappa_rho,
                  bar_epsilon_1_std=bar_epsilon_1_std,
                  bar_epsilon_1_mean=bar_epsilon_1_mean,
                  alpha_pi=alpha_pi)
    
    monitor=smc_monitor_f( gs=gs,
                fixed_tau=fixed_tau,
                fixed_bar_epsilon=fixed_bar_epsilon,
                constrained_epsilon_matrix=constrained_epsilon_matrix,
                relax_rho=relax_rho)
    
    if(is.null(inits)){
      inits=smc_inits(v=dim(n_vsa)[1],
               s=dim(n_vsa)[1],
               g=G,
               i=if(mcmc){1}else{n_chains},
               tau_vgb_0=tau_vgb_0,
               tau_ivgb_0=tau_ivgb_0,#set initial value
               pi_gs_0=pi_gs_0,
               pi_igs_0=pi_igs_0,#set initial value
               gs=gs,
               epsilon_ba_0=epsilon_ba_0,
               fixed_pi=fixed_pi,
               fixed_tau=fixed_tau,
               alpha_tau=alpha_tau,
               alpha_epsilon=alpha_epsilon,
               alpha_bar_epsilon=alpha_bar_epsilon,
               fixed_bar_epsilon=fixed_bar_epsilon,
               kappa_rho=kappa_rho,
               alpha_rho=alpha_rho,
               bar_epsilon_1_0=bar_epsilon_1_0,
               bar_epsilon_1_std=bar_epsilon_1_std,
               bar_epsilon_1_mean=bar_epsilon_1_mean,
               alpha_pi=alpha_pi) }  

    
  # Compiling and producing posterior samples from the model.
    
    
    
    
    
    if(gs=="custom"){
      smc_samples<-
        smc_custom(n_vsa=n_vsa,
                         tau_vgb=tau_vgb,
                         tau_ivgb_0=tau_ivgb_0,
                         tau_vgb_0=tau_vgb_0,
                         G=G,
                         pi_igs_0=pi_igs_0,#set initial value
                         pi_gs_0=pi_gs_0,#set initial value
                         block_tau=block_tau,
                         alpha_tau=alpha_tau,
                         alpha_rho=alpha_rho,
                         alpha_epsilon=alpha_epsilon,
                         bar_epsilon_1_std=bar_epsilon_1_std,
                         bar_epsilon_1_mean=bar_epsilon_1_mean,
                         alpha_bar_epsilon=alpha_bar_epsilon,
                         bar_epsilon_1=bar_epsilon_1,
                         kappa_rho=kappa_rho,
                         alpha_pi=alpha_pi,
                         n_chains = n_chains,
                         n_vsa_df=n_vsa_df,
                         g=G,
                         t_min=t_min,
                         t_max=t_max,
                         ess_min=ess_min,
                         trace_all=trace_all,
                         .update_lambda=update_lambda,
                         mcmc=mcmc,
                         inits=inits,
                         tempering_n=tempering_n,
                         bar_epsilon_1_tempering=bar_epsilon_1_tempering,
                         alpha_tau_tempering=alpha_tau_tempering,
                         alpha_tau_seq=alpha_tau_seq,
                         bar_epsilon_1_seq=bar_epsilon_1_seq)
      }
    
   if(gs=="jags"){
     smc_samples <-runjags::run.jags(
      model = model_string,
      data = observations_and_constants,
      monitor = monitor,n.chains = n_chains,
      #inits=inits,
      ...
     )}
    if(gs=="nimble"){
      nimble_code=eval(parse(text=paste0("nimble::nimbleCode(",model_string,")")))
      model <- nimble::nimbleModel(
        code = nimble_code, 
        data = observations_and_constants$data_list, 
        constants = observations_and_constants$constants, 
        inits = inits)
      compiled_model <- compileNimble(model)   
      mcmcConf <- configureMCMC(model)
      Rmcmc <- buildMCMC(mcmcConf)
      Cmodel <- compileNimble(model)
      Cmcmc <- compileNimble(Rmcmc, project = model)
      
      smc_samples <- runMCMC(Cmcmc, ...)
      
    }
    
    
    
    
    if(gs=="stan"){
      smc_samples=rstan::stan(model_code = model_string, 
                              data = observations_and_constants, ...)}
    
    
   list(smc_samples=smc_samples,
        model_string=model_string,
        monitor=monitor)
  
 }
  
