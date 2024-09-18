#' Sequential Monte Carlo (SMC) or MCMC Sampling for Bayesian Inference
#'
#' This function runs a Sequential Monte Carlo (SMC) or Markov Chain Monte Carlo (MCMC)
#' algorithm to sample from a model specified in JAgS, NIMBLE, Stan, or a custom model.
#' It supports various model structures with parameters for block sampling of tau, 
#' constrained epsilon matrix, and tempering of parameters. 
#'
#' @param n_vsa Matrix or array of observed data.
#' @param gs String, the model specification engine. Options are "jags", "nimble", "stan", or "custom".
#' @param tau_ivgb_0 Initial values for tau_ivgb (optional).
#' @param tau_vgb_0 Initial values for tau_vgb (optional).
#' @param pi_igs_0 Initial values for pi_igs (optional).
#' @param pi_gs_0 Initial values for pi_gs (optional).
#' @param tau_vgb Matrix, initial values for tau_vgb. If NULL, it will be estimated.
#' @param g Integer, number of variants (automatically determined if tau_vgb is provided).
#' @param block_tau Logical, whether to use block sampling for tau. Defaults to TRUE.
#' @param alpha_tau Numeric, prior parameter for tau (optional).
#' @param alpha_epsilon Numeric, prior parameter for epsilon (optional).
#' @param bar_epsilon_1_std Numeric, standard deviation for bar_epsilon_1 (optional).
#' @param bar_epsilon_1_mean Numeric, mean for bar_epsilon_1 (optional).
#' @param alpha_bar_epsilon Numeric vector, prior parameters for bar_epsilon (default: c(1, 10)).
#' @param bar_epsilon_1 Initial value for bar_epsilon_1 (optional).
#' @param kappa_rho Numeric, prior parameter for rho (optional).
#' @param alpha_rho Numeric, prior parameter for rho (optional).
#' @param alpha_pi Numeric, prior parameter for pi (default: 0.1).
#' @param n_chains Integer, number of chains for MCMC sampling (default: 2).
#' @param mcmc Logical, whether to use MCMC (default: FALSE for SMC).
#' @param inits List, initial values for the parameters (optional).
#' @param t_min Integer, minimum temperature for tempering (default: 1).
#' @param t_max Integer, maximum temperature for tempering (default: 30).
#' @param ess_min Minimum effective sample size (optional).
#' @param smc_kernel Function, SMC kernel to use (default: `desman_kernel`).
#' @param trace_all Logical, whether to trace all samples (default: TRUE).
#' @param .update_lambda Function, function to update lambda (optional).
#' @param n_vsa_df Data frame, additional data for n_vsa (optional).
#' @param tempering_n Logical, whether to apply tempering to the number of samples (default: FALSE).
#' @param bar_epsilon_1_tempering Logical, whether to apply tempering to bar_epsilon_1 (default: FALSE).
#' @param alpha_tau_tempering Logical, whether to apply tempering to alpha_tau (default: FALSE).
#' @param alpha_tau_seq Numeric vector, sequence for tempering alpha_tau (optional).
#' @param bar_epsilon_1_seq Numeric vector, sequence for tempering bar_epsilon_1 (optional).
#' @param ... Additional arguments to pass to the model sampler (JAgS, Stan, NIMBLE, or custom).
#' 
#' @return A list containing the SMC or MCMC samples (`smc_samples`), the model string (`model_string`), and the monitor (`monitor`).
#' 
#' @export
#'@examples
#'gs="custom"
#'tau_pi_n <- sim_tau_pi_epsilon_n(v = 50, g = 5, s = 3, n = 1000, alpha_pi = 1)
#'n_vsa = tau_pi_n$n_vsa
#'tau_vgb = tau_pi_n$tau_vgb
#'g=if(!is.null(tau_vgb)){dim(tau_vgb)[2]}else{5}
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
#'g=g,
#'gs=gs,
#'adapt=0,
#'burnin=0,
#'sample=10)
#')}
#'X|>mcmc_output_df("bar_epsilon")|>View()
    
smc_run<-
  function(n_vsa,
           gs="custom",
           tau_ivgb_0=NULL,
           tau_vgb_0=NULL,
           pi_igs_0=NULL,#set initial value
           pi_gs_0=NULL,#set initial value
           tau_vgb=NULL,
           gstar=if(!is.null(tau_vgb)){dim(tau_vgb)[2]}else{NULL},
           g=5,
           block_tau=FALSE,
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
    print("smc_run1")
    print(g)
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
                  g=g,
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
               s=dim(n_vsa)[2],
               g=g,
               i=n_chains,
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
                         g=g,
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
  
