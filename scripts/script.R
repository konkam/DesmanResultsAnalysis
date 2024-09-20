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
set.seed(1)
tau_pi_n <- sim_tau_pi_epsilon_n(v = 4, g = gm, s = 3, 
n = 1000, alpha_pi = .1, bar_epsilon_1 = .001)
n_vsa =tau_pi_n$n_vsa
smc_run(n_vsa,
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


###################
g=2;v=2;s=3;alpha_pi=2
tau_vgb=c("aa","cc")|>translate_dna_string_vector_to_string_matrix()|>translate_dna_matrix_to_binary_array()
tau_vgb0=c("ac","ca")|>translate_dna_string_vector_to_string_matrix()|>translate_dna_matrix_to_binary_array()

pi_gs <- sim_pi_gs(g = g, s = s, alpha_pi = alpha_pi)
reorder_g=order(pi_gs|>plyr::aaply(1,sum),decreasing = TRUE)|>
  (function(x){x|>setNames(x)})()
pi_gs=pi_gs[reorder_g,]|>(`dimnames<-`)(list(g=1:g,s=1:s))
tau_vgb=tau_vgb[,reorder_g,]|>(`dimnames<-`)(list(v=1:v,g=1:g,b=nucleotides))
epsilon_ba=epsilon_ba_f(.01)
n_vsa <- sim_n_vsa(n = 1000, tau_vgb = tau_vgb, pi_gs = pi_gs, epsilon_ba = epsilon_ba)
inits=plyr::rlply(20,list(tau_vgb=tau_vgb0,pi_gs =sim_pi_gs(g = g, s = s, alpha_pi = alpha_pi)))
smc_run(n_vsa,
       gs="jags",
       tau_vgb=NULL,
       G=g,
       block_tau=FALSE,
       alpha_tau=NULL,
       alpha_epsilon=NULL,
       bar_epsilon_1_std=NULL,
       bar_epsilon_1_mean=NULL,
       alpha_bar_epsilon=NULL,
       bar_epsilon_1=0.001,
       kappa_rho=NULL,
       alpha_rho=NULL,
       alpha_pi=.1,
       n_chains = 20,
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
true_parameter2=list(tau_vgb=tau_vgb,pi_gs=pi_gs)
smc_samples2|>
  plot_pi_tau_from_sample(n_vsa,v=1:3,
                          importance_g_it=importance_g_it2,
                          reorder_g=reorder_g2,
                          true_parameter=true_parameter2,t_step=1)

###################



set.seed(1)
tau_pi_n <- sim_tau_pi_epsilon_n(v = 4, g = gm, s = 3, 
                                 n = 10, alpha_pi = .1, bar_epsilon_1 = .001)
n_vsa =tau_pi_n$n_vsa


smc_run(n_vsa,
       gs="jags",
       tau_vgb=NULL,
       G=5,
       block_tau=FALSE,
       alpha_tau=.001,
       alpha_epsilon=NULL,
       bar_epsilon_1_std=NULL,
       bar_epsilon_1_mean=NULL,
       alpha_bar_epsilon=NULL,
       bar_epsilon_1=.2,
       kappa_rho=NULL,
       alpha_rho=NULL,
       alpha_pi=1,
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
