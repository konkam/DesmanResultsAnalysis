library(tidyverse)
library(targets)
library(tarchetypes)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(runjags)
library(RColorBrewer)

"R" |>  list.files(full.names = TRUE) |>sapply(FUN = source)

v=2;s=2;G=g=2
n_chains=i=3
n_vsa=abind::abind(
  matrix(c(600,400,0,1,400,600,1,0),2,4,byrow=TRUE),
  matrix(c(800,200,0,1,200,800,0,1),2,4,byrow=TRUE),along=3)|>
  aperm(c(1,3:2))|>
  namedims(index="vsa")

pi_igs_0=plyr::raply(i,matrix(c(.6,.4,.8,.2),2,2))
epsilon_iba=epsilon_iba_f(3,0.001)


tau_ivgb_0=plyr::llply(list(c("aa","cc"),c("ac","ca"),c("gg","tt")),translate_dna_string_vector_to_string_matrix)|>
  plyr::laply(translate_dna_matrix_to_binary_array)|>(function(x){names(dimnames(x))[1]="i";dimnames(x)[1]=list(1:3);x})()

desman_kernel=kernel_f(fixed_bar_epsilon=FALSE,
         constrained_epsilon_matrix=FALSE,
         block_tau=FALSE,
         fixed_tau=FALSE,
         relax_tau=FALSE,
         relax_rho=FALSE)


AM_1kernel=kernel_f(fixed_bar_epsilon=FALSE,
                           constrained_epsilon_matrix=FALSE,
                           block_tau=FALSE,
                           fixed_tau=FALSE,
                           relax_tau=FALSE,
                           relax_rho=FALSE)


inits=smc_inits(i = i,s=s,v=v,g=g,tau_ivgb_0 =tau_ivgb_0,pi_igs_0 =pi_igs_0)
desman_inits=c(inits,
               list(epsilon_iba=epsilon_iba_f(i,.001)))


desman=smc_custom(n_vsa,
           gs="custom",
           G=2,
           pi_igs_0=pi_igs_0,#set initial value
           block_tau=FALSE,
           bar_epsilon_1 = NULL,
           alpha_epsilon=.1,
           alpha_pi=.1,
           n_chains = 3,
           t_max=30,
           mcmc=TRUE,
           inits=desman_inits)





desman$theta$tau_ivgb|>
  plyr::aaply(1,translate_dna_binary_array_to_string_vector)

desman$theta$pi_igs

AM1=smc_custom(n_vsa,
                    gs="custom",
                    G=2,
                    pi_igs_0=pi_igs_0,#set initial value
                    block_tau=TRUE,
                    alpha_tau=NULL,
                    alpha_rho=NULL,
                    alpha_epsilon=NULL,
                    bar_epsilon_1_std=.01,
                    bar_epsilon_1_mean=.001,
                    alpha_bar_epsilon=NULL,
                    bar_epsilon_1=NULL,
                    kappa_rho=NULL,
                    alpha_pi=.1,
                    n_chains = 3,
                    n_vsa_df=NULL,
                    g=G,
                    t_min=1,
                    t_max=30,
                    ess_min=NULL,
                    smc_kernel=desman_kernel,
                    trace_all=TRUE,
                    .update_lambda=update_lambda,
                    mcmc=FALSE,
                    inits=inits,
                    tempering_n=FALSE,
                    bar_epsilon_1_tempering=FALSE,
                    alpha_tau_tempering=FALSE,
                    alpha_tau_seq=NULL,
                    bar_epsilon_seq=NULL)
    
AM1$theta$tau_ivgb|>
  plyr::aaply(1,translate_dna_binary_array_to_string_vector)

AM1$theta$pi_igs


AM2=smc_custom(n_vsa,
               gs="custom",
               G=2,
               tau_vgb = c("ac","ca")|>translate_dna_string_vector_to_string_matrix()|>translate_dna_matrix_to_binary_array(),
               block_tau=FALSE,
               bar_epsilon_1_std=.01,
               bar_epsilon_1_mean=.001,
               alpha_bar_epsilon=NULL,
               bar_epsilon_1=NULL,
               kappa_rho=NULL,
               alpha_pi=.1,
               n_chains = 3,
               n_vsa_df=NULL,
               g=G,
               t_min=1,
               t_max=30,
               ess_min=NULL,
               smc_kernel=desman_kernel,
               trace_all=TRUE,
               .update_lambda=update_lambda,
               mcmc=TRUE,
               tempering_n=FALSE,
               bar_epsilon_1_tempering=FALSE,
               alpha_tau_tempering=FALSE,
               alpha_tau_seq=NULL,
               bar_epsilon_seq=NULL)


AM2$theta$pi_igs

AM3=smc_custom(n_vsa,
               gs="custom",
               G=2,
               tau_vgb = NULL,
               tau_vgb_0=tau_vgb_0,
               block_tau=FALSE,
               bar_epsilon_1_std=NULL,
               bar_epsilon_1_mean=NULL,
               alpha_bar_epsilon=NULL,
               bar_epsilon_1=NULL,
               alpha_rho=.001,
               alpha_pi=.1,
               n_chains = 3,
               n_vsa_df=NULL,
               g=G,
               t_min=1,
               t_max=30,
               ess_min=NULL,
               trace_all=TRUE,
               mcmc=TRUE,
               tempering_n=FALSE,
               bar_epsilon_1_tempering=FALSE,
               alpha_tau_tempering=FALSE,
               alpha_tau_seq=NULL,
               bar_epsilon_seq=NULL)


AM3$theta$pi_igs


AM3$theta$tau_ivgb|>
  plyr::aaply(1,translate_dna_binary_array_to_string_vector)


AM4=smc_custom(n_vsa,
               gs="custom",
               G=2,
               tau_vgb = NULL,
               tau_vgb_0=tau_vgb_0,
               block_tau=FALSE,
               bar_epsilon_1_std=.01,
               bar_epsilon_1_mean=.001,
               alpha_bar_epsilon=NULL,
               bar_epsilon_1=NULL,
               alpha_rho=NULL,
               alpha_pi=1,
               n_chains = 3,
               n_vsa_df=NULL,
               g=G,
               alpha_tau=.001,
               t_min=1,
               t_max=30,
               ess_min=NULL,
               trace_all=TRUE,
               mcmc=TRUE,
               tempering_n=FALSE,
               bar_epsilon_1_tempering=FALSE,
               alpha_tau_tempering=FALSE,
               alpha_tau_seq=NULL,
               bar_epsilon_seq=NULL)


AM4$theta$pi_igs


AM4$theta$tau_ivgb|>
  plyr::aaply(1:3,function(x){1*(x==max(x))})|>
  plyr::aaply(1,translate_dna_binary_array_to_string_vector)




AM4_temp=smc_custom(n_vsa,
               gs="custom",
               G=2,
               tau_vgb = NULL,
               tau_vgb_0=tau_vgb_0,
               block_tau=FALSE,
               bar_epsilon_1_std=.01,
               bar_epsilon_1_mean=.001,
               alpha_bar_epsilon=NULL,
               bar_epsilon_1=NULL,
               alpha_rho=NULL,
               alpha_pi=1,
               n_chains = 3,
               n_vsa_df=NULL,
               g=G,
               alpha_tau=.001,
               t_min=1,
               t_max=t_max,
               ess_min=NULL,
               trace_all=TRUE,
               mcmc=TRUE,
               tempering_n=FALSE,
               bar_epsilon_1_tempering=FALSE,
               alpha_tau_tempering=TRUE,
               alpha_tau_seq=max(.001,0.499*(1-(1:t_max)/t_max)^(1/4)),
               bar_epsilon_seq=NULL)


AM4$theta$pi_igs


AM4$theta$tau_ivgb|>
  plyr::aaply(1:3,function(x){1*(x==max(x))})|>
  plyr::aaply(1,translate_dna_binary_array_to_string_vector)

