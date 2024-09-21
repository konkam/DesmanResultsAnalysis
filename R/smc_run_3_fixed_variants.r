#' This function performs Sequential Monte Carlo (SMC) sampling with fixed tau. 
#' It adjusts the `tau_vgb` based on the given `n_vsa` and other parameters,
#' then calls the `smc_run` function to complete the process.
#'
#' @param n_vsa Integer. The number of VSA (variable selection algorithm) samples to run.
#' @param tau_vgb Matrix or list. The tau matrix related to VGB (variance-gamma model parameters).
#' @param discriminant_v vector of positions index for which tau_vgb is not constant, default is computed from `discriminant_v_f(tau_vgb)`.
#' @param reduced_tau_vgb List. A reduced version of `tau_vgb`, computed by default using `reduce_tau_vgb_f(tau_vgb, n_vsa, discriminant_v)`.
#' @param g Integer. Number of groups or dimensions, derived from `reduced_tau_vgb$tau_vgb`. Default is the number of columns in `reduced_tau_vgb$tau_vgb`.
#' @param alpha_pi Numeric. Alpha prior for Pi, default is 0.1.
#' @param alpha_bar_epsilon Numeric vector. Prior parameters for epsilon, default is c(1, 10).
#' @param n_chains Integer. The number of Markov chains to run, default is 2.
#' @param ... Additional arguments passed to the `smc_run` function.
#'
#' @return The output from the `smc_run` function, which may include posterior samples or SMC diagnostics.
#'
#' @details
#' This function prepares and processes the `tau_vgb` variable by reducing it 
#' based on the discriminant vector and the number of VSA samples. It then runs
#' the Sequential Monte Carlo algorithm using `smc_run`, passing necessary parameters.
#' If the number of groups `g` matches the dimensions of `reduced_tau_vgb$tau_vgb`, 
#' an initial Pi vector (`pi_igs_0`) is set up.
#'
#' @examples
#'gs="custom"
#'tau_pi_n <- sim_tau_pi_epsilon_n(v = 50, g = 5, s = 3, n = 1000, alpha_pi = 1)
#'n_vsa = tau_pi_n$n_vsa
#'tau_vgb = tau_pi_n$tau_vgb
#'G=if(!is.null(tau_vgb)){dim(tau_vgb)[2]}else{5}
#'bar_epsilon=NULL
#'alpha_tau=NULL
#'alpha_bar_epsilon=c(1,10)
#'bar_epsilon_1_std=NULL
#'bar_epsilon_1_mean=NULL
#'alpha_pi=NULL
#'fixed_bar_epsilon=!is.null(bar_epsilon_1)
#'constrained_epsilon_matrix=TRUE
#'n_chains = 2
#'list(n_vsa=n_vsa,
#'tau_vgb=tau_vgb,
#'G=G,
#'alpha_pi=alpha_pi,
#'block_tau=FALSE,
#'alpha_tau=NULL,
#'alpha_epsilon=NULL,
#'bar_epsilon_1=NULL,
#'alpha_bar_epsilon=alpha_bar_epsilon,
#'kappa_rho=NULL,
#'n_chains = n_chains)|>list2env(.GlobalEnv)
#'tau_vgb=NULL
#'tar_load_everything()
#'formals(smc_run)|>
#'(function(x){
#'x[sapply(x,function(xx){!is.element("name",class(xx))})]})()|>
#'plyr::llply(eval)|>
#'list2env(.GlobalEnv)
#'list(n_vsa=n_vsa,
#'      tau_vgb=tau_vgb,
#'      G=G,
#'     alpha_pi=alpha_pi,
#'     block_tau=FALSE,
#'     alpha_tau=NULL,
#'      alpha_epsilon=NULL,
#'     bar_epsilon_1=NULL,
#'     alpha_bar_epsilon=alpha_bar_epsilon,
#'      kappa_rho=NULL,
#'      n_chains = n_chains)|>
#'list2env(.GlobalEnv)
       
       

#'smc_run_fixed_tau(n_vsa,
#'tau_vgb=tau_pi_n$tau_vgb,  
#'G=G,
#'gs=gs,
#'alpha_bar_epsilon=alpha_bar_epsilon,
#'bar_epsilon_1_std=bar_epsilon_1_std,
#'bar_epsilon_1_mean=bar_epsilon_1_mean,
#'alpha_pi=1,
#'n_chains = 2,
#'burnin=4,
#'adapt=3,
#'sample=10)->X
#'X|>mcmc_output_df(variable_name = "pi")|>View()
#'X|>mcmc_output_df(variable_name = "bar_epsilon")|>View()
#'X|>mcmc_output_df(variable_name = "pi",index = "gs")|>
#'ggplot(aes(x=iteration,y=value,group=interaction(chain,g,s),fill=g))+
#'geom_area()+
#'facet_grid(s~chain)
#' @export


smc_run_fixed_tau <- function(n_vsa,
                            tau_vgb,
                            discriminant_v=discriminant_v_f(tau_vgb),
                            reduced_tau_vgb=reduce_tau_vgb_f(tau_vgb,n_vsa,discriminant_v),
                            g=dim(reduced_tau_vgb$tau_vgb)[2],
                            alpha_pi=.1,
                            alpha_bar_epsilon=c(1,10),
                            n_chains =2,
                            ...) {
  
  
  pi_igs_0=if(g==dim(reduced_tau_vgb$tau_vgb)[2]){plyr::raply(n_chains,reduced_tau_vgb$pi_gs_0)}else{NULL}
  
  print("smc_run_fixed_tau")
  smc_run(n_vsa=n_vsa,
         tau_vgb=reduced_tau_vgb$tau_vgb,
         g=g,
         alpha_pi=alpha_pi,
         alpha_bar_epsilon=alpha_bar_epsilon,
         n_chains = n_chains,
         pi_igs_0 = pi_igs_0,
         ...) 
  }
