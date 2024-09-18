#'@examples
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


smc_run_fixed_tau <- function(n_vsa,
                            tau_vgb,
                            g=dim(tau_vgb)[],
                            alpha_pi=.1,
                            alpha_bar_epsilon=c(1,10),
                            n_chains =2,
                            ...) {
  smc_run(n_vsa=n_vsa,
         tau_vgb=tau_vgb,
         g=g,
         alpha_pi=alpha_pi,
         alpha_bar_epsilon=alpha_bar_epsilon,
         n_chains = n_chains,
         ...) 
  }
