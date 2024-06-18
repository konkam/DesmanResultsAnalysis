#'@examples
#'gs="jags"
#'tau_pi_n <- sim_tau_pi_epsilon_n(v = 3, g = 3, s = 5, 
#'n = 10000, alpha_pi = 1, epsilon_bar_1 = .001)
#'n_vsa = tau_pi_n$n_vsa
#'G=3
#'
#'inference_5_run(n_vsa,
#'G=G,
#'gs=gs,
#'alpha_tau=.2,
#'bar_epsilon_1_std=.1,
#'bar_epsilon_1_mean=.01,
#'alpha_pi=1,
#'n_chains = 4,
#'burnin=0,
#'adapt=0,
#'sample=1000,
#'thin=1)->smc_sample
#'smc_sample$smc_samples->smc_samples
#'reorder_g=smc_reorder(smc_samples,n_vsa)
#'importance_g_it=
#'  importance_g_it_f(smc_samples=smc_samples,
#'                    n_vsa=n_vsa,
#'                    reorder_g=reorder_g)
#'true_parameter=tau_pi_n
#'smc_samples|>
#'plot_pi_tau_from_sample(n_vsa,v=1:3,
#'importance_g_it=importance_g_it,
#'reorder_g=reorder_g,
#'true_parameter=tau_pi_n,t_step=1)
#'tau_pi_n$tau_vgb[1:3,,]|>translate_dna_binary_array_to_string_vector()
#'tau_pi_n$pi_gs|>plyr::aaply(1,sum)
#'tau_pi_n$pi_gs
#'tau_pi_n$n_vsa[1:2,,]|>plyr::aaply(c(1,3),sum)

inference_5_run <- function(n_vsa,
                            G,
                            gs="jags",
                            alpha_pi=.1,
                            alpha_tau=.001,
                            bar_epsilon_1_std=NULL,
                            bar_epsilon_1_mean=NULL,
                            alpha_bar_epsilon=c(1,100),
                            n_chains =2,
                            ...) {
  gs_run(n_vsa,
         gs=gs,
         tau_vgb=NULL,
         G=G,
         block_tau=FALSE,
         alpha_tau=alpha_tau,
         alpha_epsilon=NULL,
         bar_epsilon_1_std=bar_epsilon_1_std,
         bar_epsilon_1_mean=bar_epsilon_1_mean,
         alpha_bar_epsilon=alpha_bar_epsilon,
         bar_epsilon=NULL,
         kappa_rho=NULL,
         alpha_pi=alpha_pi,
         n_chains = n_chains,
         ...) }
