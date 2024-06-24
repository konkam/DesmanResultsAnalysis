#'@examples
#'gs="jags"
#'tau_pi_n <- sim_tau_pi_epsilon_n(v = 50, g = 5, s = 3, 
#'n = 1000, alpha_pi = .1, bar_epsilon_1 = .000001)
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
#'
#'inference_4_run(n_vsa,
#'G=G,
#'gs=gs,
#'kappa_rho=c(1,1000),
#'alpha_pi=.1,
#'n_chains = 2,
#'burnin=1004,
#'adapt=3000,
#'sample=1000)->X
#'X|>mcmc_output_df(variable_name = "rho",index="vga")|>
#'dplyr::filter(v==1)|>
#'ggplot(aes(x=iteration,y=value,group=interaction(chain,g,a),fill=a))+
#'geom_area()+
#'facet_grid(g~chain)
#'tau_pi_n$tau_vgb[1,,]  
#'
#'X|>mcmc_output_df(variable_name = "pi",index = "gs")|>
#'ggplot(aes(x=iteration,y=value,group=interaction(chain,g,s),fill=g))+
#'geom_area()+
#'facet_grid(s~chain)
#'tau_pi_n$pi_gs  
#'

inference_4_run <- function(n_vsa,
                            G,
                            tau_vgb,
                            gs="jags",
                            alpha_pi=.1,
                            kappa_rho=c(1,100),
                            n_chains =2,
                            ...) {
  gs_run(n_vsa,
         gs=gs,
         tau_vgb=NULL,
         G=G,
         block_tau=FALSE,
         alpha_tau=NULL,
         alpha_epsilon=NULL,
         bar_epsilon_1_std=NULL,
         bar_epsilon_1_mean=NULL,
         alpha_bar_epsilon=NULL,
         bar_epsilon=NULL,
         kappa_rho=kappa_rho,
         alpha_pi=alpha_pi,
         n_chains = n_chains,
         ...) }
