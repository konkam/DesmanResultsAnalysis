library(tidyverse)
library(targets)
library(tarchetypes)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(runjags)
library(RColorBrewer)

"R" |>  list.files(full.names = TRUE) |>sapply(FUN = source)
gm=4
g=4


gs="jags"
tau_pi_n <- sim_tau_pi_epsilon_n(v = 4, g = gm, s = 3, 
n = 1000, alpha_pi = .1, epsilon_bar_1 = .001)
n_vsa =tau_pi_n$n_vsa
gs_run(n_vsa,
       gs="jags",
       tau_vgb=NULL,
       G=g,
       block_tau=TRUE,
       alpha_tau=NULL,
       alpha_epsilon=NULL,
       bar_epsilon_1_std=NULL,
       bar_epsilon_1_mean=NULL,
       alpha_bar_epsilon=NULL,
       bar_epsilon_1=.001,
       kappa_rho=NULL,
       alpha_rho=NULL,
       alpha_pi=.1,
       n_chains = 4,
       burnin=0,
       adapt=0,
       sample=40,
       thin=500)->smc_sample
       smc_sample$smc_samples->smc_samples
reorder_g=smc_reorder(smc_samples,n_vsa)
importance_g_it=
  importance_g_it_f(smc_samples=smc_samples,
                    n_vsa=n_vsa,
                    reorder_g=reorder_g)
true_parameter=tau_pi_n
smc_samples|>
plot_pi_tau_from_sample(n_vsa,v=1:3,
importance_g_it=importance_g_it,
reorder_g=reorder_g,
true_parameter=tau_pi_n,t_step=1)

tau_pi_n$tau_vgb|>translate_dna_binary_array_to_string_vector()
tau_pi_n$pi_gs|>plyr::aaply(1,sum)
tau_pi_n$pi_gs
n_vsa




gs_run(n_vsa,
       gs="jags",
       tau_vgb=NULL,
       G=g,
       block_tau=FALSE,
       alpha_tau=NULL,
       alpha_epsilon=NULL,
       bar_epsilon_1_std=NULL,
       bar_epsilon_1_mean=NULL,
       alpha_bar_epsilon=NULL,
       bar_epsilon_1=.001,
       kappa_rho=NULL,
       alpha_rho=NULL,
       alpha_pi=.1,
       n_chains = 4,
       burnin=0,
       adapt=0,
       sample=40,
       thin=500)->smc_sample2
smc_sample2$smc_samples->smc_samples2
reorder_g2=smc_reorder(smc_samples2,n_vsa)
importance_g_it2=
  importance_g_it_f(smc_samples=smc_samples2,
                    n_vsa=n_vsa,
                    reorder_g=reorder_g2)
true_parameter=tau_pi_n
smc_samples2|>
  plot_pi_tau_from_sample(n_vsa,v=1:3,
                          importance_g_it=importance_g_it2,
                          reorder_g=reorder_g2,
                          true_parameter=tau_pi_n,t_step=1)

###################




gs_run(n_vsa,
       gs="jags",
       tau_vgb=NULL,
       G=12,
       block_tau=FALSE,
       alpha_tau=NULL,
       alpha_epsilon=NULL,
       bar_epsilon_1_std=NULL,
       bar_epsilon_1_mean=NULL,
       alpha_bar_epsilon=NULL,
       bar_epsilon_1=.2,
       kappa_rho=NULL,
       alpha_rho=NULL,
       alpha_pi=.1,
       n_chains = 4,
       burnin=0,
       adapt=0,
       sample=40,
       thin=500)->smc_sample3
smc_sample3$smc_samples->smc_samples3
reorder_g3=smc_reorder(smc_samples3,n_vsa)
importance_g_it3=
  importance_g_it_f(smc_samples=smc_samples3,
                    n_vsa=n_vsa,
                    reorder_g=reorder_g3)
true_parameter=tau_pi_n
smc_samples3|>
  plot_pi_tau_from_sample(n_vsa,v=1:3,
                          importance_g_it=importance_g_it3,
                          reorder_g=reorder_g3,
                          true_parameter=tau_pi_n,t_step=1)
