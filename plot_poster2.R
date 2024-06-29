library(tidyverse)
library(targets)
library(tarchetypes)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(runjags)
library(RColorBrewer)

"R" |>  list.files(full.names = TRUE) |>sapply(FUN = source)
t_max=300
v=2;s=2;G=g=2
n_chains=i=3
n_vsa=abind::abind(
  matrix(c(500,500,0,1,500,500,1,0),2,4,byrow=TRUE),
  matrix(c(500,500,0,1,500,500,0,1),2,4,byrow=TRUE),along=3)|>
  aperm(c(1,3:2))|>
  namedims(index="vsa")
tau_vgb=array(c(1,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0),c(2,2,4))|>namedims("vgb")
pi_gs=matrix(c(.5,.5,.5,.5),2,2)|>namedims("gs")

pi_igs_0=plyr::raply(i,pi_gs)|>namedims("igs")
epsilon_ba=epsilon_ba_f(0.001)|>namedims("ba")
epsilon_iba=epsilon_iba_f(3,0.001)
true_parameter=list(tau_vgb=tau_vgb,epsilon_ba=epsilon_ba,pi_gs=pi_gs)

tau_ivgb_0=plyr::llply(list(c("aa","cc"),c("ac","ca"),c("gg","tt")),translate_dna_string_vector_to_string_matrix)|>
  plyr::laply(translate_dna_matrix_to_binary_array)|>(function(x){names(dimnames(x))[1]="i";dimnames(x)[1]=list(1:3);x})()



inits=smc_inits(i = i,s=s,v=v,g=g,tau_ivgb_0 =tau_ivgb_0,pi_igs_0 =pi_igs_0)
desman_inits=c(inits,
               list(epsilon_iba=epsilon_iba_f(i,.001)))

set.seed(3)
desman=smc_custom(n_vsa,
                  gs="custom",
                  G=2,
                  tau_ivgb_0=tau_ivgb_0,
                  pi_igs_0=pi_igs_0,#set initial value
                  block_tau=FALSE,
#                  bar_epsilon_1 = NULL,
                  bar_epsilon_1=.001,
#                  alpha_epsilon=.1,
                  alpha_pi=1,
                  n_chains = 3,
                  t_max=t_max,
                  mcmc=TRUE,
                  inits=desman_inits)
plot_desman=plot_pi_tau_from_sample(
  smc_samples=NULL,
  traces=desman$traces,
  what=c("tau"),
  n_vsa=n_vsa,
  t_step = 1,
  true_parameter = NULL,
  reorder_g=data.frame(g=1:g,reorder_g=1:g))$plot_tau

ggsave(plot=plot_desman,path ="latex/poster/",
       filename= "plotdesman.pdf",height = 10,width=15,units="cm")

trueposterior_tau=list(tau_tivgb=
    plyr::laply(list(c("ac","ca"),c("aa","cc")),
                            translate_dna_string_vector_to_string_matrix)|>
         plyr::aaply(1,translate_dna_matrix_to_binary_array)|>
      array(c(1,2,2,2,4))|>
              namedims('tivgb'),
    pi_tigs=
      array(.5,c(1,2,2,2))|>
      namedims('tigs'))

plotpt=convert_custom_ouput(trueposterior_tau)$tau_vgb|>
  (function(x){rbind(x,x|>dplyr::mutate(iteration=2))})()|>
  dplyr::mutate(g=paste0("g=",g),v=paste0("v=",v))|>
  ggplot(aes(x=iteration,y=value,group=interaction(chain,g,b,v),fill=b))+
  geom_area()+
  xlab("")+ylab("")+
    facet_grid(g+v~chain,scales="free_x")+ 
  theme(axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()) 
  ggsave(plot=plotpt,path ="latex/poster/",
         filename= "posteriortrue.pdf",height = 10,width=15,units="cm")


desman$theta$tau_ivgb|>
  plyr::aaply(1,translate_dna_binary_array_to_string_vector)


desman$theta$tau_ivgb|>
  plyr::aaply(1,translate_dna_binary_array_to_string_vector)

desman$theta$pi_igs
set.seed(1)
AM1_inits=c(inits,list(epsilon_iba=epsilon_iba_f(i=3,bar_epsilon_1=,.001)))
AM1=smc_custom(n_vsa,
               gs="custom",
               G=2,
               tau_vga=NULL,
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
               alpha_pi=1,
               n_chains = 3,
               n_vsa_df=NULL,
               g=G,
               t_min=1,
               t_max=t_max,
               ess_min=NULL,
               smc_kernel=desman_kernel,
               trace_all=TRUE,
               .update_lambda=update_lambda,
               mcmc=FALSE,
               inits=AM1_inits,
               tempering_n=FALSE,
               bar_epsilon_1_tempering=FALSE,
               alpha_tau_tempering=FALSE,
               alpha_tau_seq=NULL,
               bar_epsilon_seq=NULL)


plotblock=plot_pi_tau_from_sample(
  smc_samples=NULL,
  traces=AM1$traces,
  n_vsa=n_vsa,
  t_step = 1,
  true_parameter = true_parameter,
  reorder_g=data.frame(g=1:g,reorder_g=1:g))$plot_tau

ggsave(plot=plotblock,path ="latex/poster/",
       filename= "plotblock.pdf",height = 10,width=15,units="cm")


set.seed(1)
AM1temp=smc_custom(n_vsa,
               gs="custom",
               G=2,
               tau_vga=NULL,
               tau_ivgb_0=NULL,
               pi_gs_0=NULL,
               pi_igs_0=NULL,#set initial value
               block_tau=FALSE,
               alpha_tau=NULL,
               alpha_rho=NULL,
               alpha_epsilon=NULL,
               bar_epsilon_1_std=NULL,
               bar_epsilon_1_mean=NULL,
               alpha_bar_epsilon=NULL,
               bar_epsilon_1=NULL,
               kappa_rho=NULL,
               alpha_pi=.1,
               n_chains = 50,
               n_vsa_df=NULL,
               g=G,
               t_min=1,
               t_max=t_max,
               ess_min=NULL,
               smc_kernel=desman_kernel,
               trace_all=TRUE,
               .update_lambda=update_lambda,
               mcmc=TRUE,
               inits=NULL,
               tempering_n=FALSE,
               bar_epsilon_1_tempering=TRUE,
               alpha_tau_tempering=FALSE,
               alpha_tau_seq=NULL,
               bar_epsilon_1_seq=(function(x){.001+0.249*(1/x)^(.1)})(1:1000))


plotepsilontemp=plot_pi_tau_from_sample(
  smc_samples=NULL,
  sel_i=c(4:6,8),
  traces=AM1temp$traces,
  n_vsa=n_vsa,
  t_step = 1,
  true_parameter = NULL,
  reorder_g=data.frame(g=1:g,reorder_g=1:g))$plot_tau

ggsave(plot=plotepsilontemp,path ="latex/poster/",
       filename= "plotepsilontemp.pdf",height = 10,width=15,units="cm")

set.seed(1)
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
               t_max=t_max,
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

plot_pi_tau_from_sample(
  smc_samples=NULL,
  traces=AM2$traces,
  what="pi",
  n_vsa=n_vsa,
  t_step = 1,
  true_parameter = true_parameter,
  reorder_g=data.frame(g=1:g,reorder_g=1:g))

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
               t_max=t_max,
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
               bar_epsilon_1_std=NULL,
               bar_epsilon_1_mean=NULL,
               alpha_bar_epsilon=NULL,
               bar_epsilon_1=0.001,
               alpha_rho=NULL,
               alpha_pi=1,
               n_chains = 3,
               n_vsa_df=NULL,
               g=G,
               alpha_tau=.01,
               t_min=1,
               t_max=t_max,
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


plot_pi_tau_from_sample(
  smc_samples=NULL,
  traces=AM4$traces,
  n_vsa=n_vsa,
  what="tau",
  t_step = 1,
  true_parameter = true_parameter,
  reorder_g=data.frame(g=1:g,reorder_g=1:g))



set.seed(1);AM4temp=smc_custom(n_vsa,
                    gs="custom",
                    G=2,
                    tau_vgb = NULL,
                    tau_vgb_0=NULL,
                    block_tau=FALSE,
                    bar_epsilon_1_std=NULL,
                    bar_epsilon_1_mean=NULL,
                    alpha_bar_epsilon=NULL,
                    bar_epsilon_1=.001,
                    alpha_rho=NULL,
                    alpha_pi=2,
                    n_chains = 30,
                    n_vsa_df=NULL,
                    g=G,
                    alpha_tau=NULL,
                    t_min=1,
                    t_max=1200,
                    ess_min=NULL,
                    trace_all=TRUE,
                    mcmc=TRUE,
                    tempering_n=FALSE,
                    bar_epsilon_1_tempering=FALSE,
                    alpha_tau_tempering=TRUE,
                    alpha_tau_seq=1-exp(-((1200:1)/1200)^4),
                    bar_epsilon_seq=NULL)


plotalpha=plot_pi_tau_from_sample(
  smc_samples=NULL,
  traces=AM4temp$traces,
  sel_i=11:20,
  n_vsa=n_vsa,
  what="tau",
  t_step = 1,
  true_parameter = true_parameter,
  reorder_g=data.frame(g=1:g,reorder_g=1:g))$plot_tau
plotalpha
ggsave(plot=plotalpha,path ="latex/poster/",
       filename= "plotalpha.pdf",height = 10,width=15,units="cm")



set.seed(1);AM4tempsmc=smc_custom(n_vsa,
                               gs="custom",
                               G=2,
                               tau_vgb = NULL,
                               tau_vgb_0=NULL,
                               block_tau=FALSE,
                               bar_epsilon_1_std=NULL,
                               bar_epsilon_1_mean=NULL,
                               alpha_bar_epsilon=NULL,
                               bar_epsilon_1=.001,
                               alpha_rho=NULL,
                               alpha_pi=2,
                               n_chains = 300,
                               n_vsa_df=NULL,
                               g=G,
                               alpha_tau=NULL,
                               t_min=1,
                               t_max=1200,
                               ess_min=NULL,
                               trace_all=TRUE,
                               mcmc=FALSE,
                               tempering_n=FALSE,
                               bar_epsilon_1_tempering=FALSE,
                               alpha_tau_tempering=TRUE,
                               alpha_tau_seq=1-exp(-((1200:1)/1200)^4),
                               bar_epsilon_seq=NULL)


plotalphasmc=plot_pi_tau_from_sample(
  smc_samples=NULL,
  traces=AM4temp$traces,
  sel_i=1:3,
  n_vsa=n_vsa,
  what="tau",
  t_step = 1,
  true_parameter = true_parameter,
  reorder_g=data.frame(g=1:g,reorder_g=1:g))$plot_tau
plotalphasmc
ggsave(plot=plotalphasmc,path ="latex/poster/",
       filename= "plotalphasmc.pdf",height = 10,width=15,units="cm")
