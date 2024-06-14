#' @description
#' Creates model string for jags or stan
#' @param bar_epsilon a numerical value. if NA, then the model uses a Dirichlet prior.
#' @param gs a character string. If "jags", the gibbs sampler used will be jags, stan if "stan".
#' @examples
#' cat(model_string_desman_f(gs="jags"))
model_string_desman_f <- function(gs="jags") {
  model_string_f(gs=gs,
                 fixed_bar_epsilon=FALSE,
                 constrained_epsilon_matrix=FALSE,
                 block_tau=FALSE,
                 fixed_tau=FALSE,
                 relax_tau=FALSE,
                 relax_rho=FALSE)}


smc_observation_desman_f<-function(n_vsa,
                          gs=gs,
                          G=G,
                          alpha_epsilon=.1,
                          alpha_pi=.1){
  
  smc_observation_generic_f(n_vsa=n_vsa,
                            gs=gs,
                            G=G,
                            alpha_bar_epsilon=NULL,
                            alpha_epsilon=alpha_epsilon,
                            bar_epsilon=NULL,
                            tau_vgb=NULL,
                            alpha_tau=NULL,
                            kappa_rho=NULL,
                            bar_epsilon_1_std=NULL,
                            bar_epsilon_1_mean=NULL,
                            alpha_pi=alpha_pi)}

#' @description
#' Run jags.
#' @param n_vsa an array of counts.
#' @param tau_vgb a collection of variants
#' @param G an integer. If G is smaller than the number of variants in the variant bin, then the algorithm is ran on the minimum of G and the number of variants in the bin.
#' @param bar_epsilon a numerical value. if NA, then the model uses a dirichlet prior.
#' @param error_rate = 0.001 controls the dirichlet prior on bar_epsilon
#' @param prior_std = 0.01 controls the dirichlet prior on bar_epsilon
#' @examples
#' tau_pi_n <- sim_tau_pi_epsilon_n(v = 50, g = 5, s = 3, n = 1000, alpha_pi = 1)
#' gs="jags"
#' block_tau=FALSE
#' X[10000,paste0("pi_gs[",c(outer(1:12,1:3,paste,sep=",")),"]")]
#' A=matrix(NA,12,3) ;for (i in 1:12){for(j in 1:3){A[i,j]=X[10000,paste0("pi_gs[",i,",",j,"]")]}}
desman_run <- function(n_vsa,
                                G,
                                gs="jags",
                                alpha_pi=.1,
                                alpha_epsilon=.1,
                                n_chains =2,
                                ...) {
  gs_run(n_vsa,
             gs="jags",
             G=G,
             tau_vgb=NULL,
             block_tau=FALSE,
             alpha_tau=NULL,
             alpha_epsilon=alpha_epsilon,
             bar_epsilon_1_std=NULL,
             bar_epsilon_1_mean=NULL,
             alpha_bar_epsilon=NULL,
             bar_epsilon=NULL,
             kappa_rho=NULL,
             alpha_pi=alpha_pi,
             n_chains = n_chains,
             ...)
  
  # Compiling and producing posterior samples from the model.
  if(gs=="jags"){
    gibbs_samples <- 
    runjags::run.jags(
      model = model_string,
      data = c(data_list,constants),
      monitor = monitor,
      ...)
    }
  
}


mcmc_desman_run <- function(n_vsa,
                            G,
                            gs="jags",
                            bar_epsilon=.01,
                            n_chains =2,
                            ...) {
  mcmc_unrelaxed_run(n_vsa,
                                                  G=G,
                                                  gs="jags",
                                                  block_tau=FALSE,
                                                  bar_epsilon_1=bar_epsilon_1,
                                                  shape_epsilon=NULL,
                                                  alpha_pi=1,
                                                  n_chains =n_chains,
                                                  ...)}


mcmc_desman_run<-function(n_vsa,
                      seed=1,
                      i,
                      g,
                      alpha_pi=.1,
                      bar_epsilon=.01,
                      shape_epsilon=c(1,100),
                      init=NULL){
  smc_sampler(n_vsa,
                        seed=1,
                        n_plus=sum(n_vsa),
                        n_vsa_df=reorder_counts(n_vsa = n_vsa,seed=seed),
                        g,
                        i,
                        t_min=1,
                        t_max,
                        ess_min,
                        smc_kernel=desman_kernel,
                        bar_epsilon=.1,
                        shape_epsilon=c(1,1000),
                        alpha_pi=.1,
                        trace_all=TRUE,
                        .update_lambda=update_lambda,
                        init=NULL,
                        mcmc=FALSE)}
