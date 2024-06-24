v=2;s=2;G=g=2
n_chains=i=3
n_vsa=abind::abind(
  matrix(c(600,400,0,1,400,600,1,0),2,4,byrow=TRUE),
  matrix(c(800,200,0,1,200,800,0,1),2,4,byrow=TRUE),along=3)|>
  aperm(c(1,3:2))|>
  namedims(index="vsa")

pi_gs_0=matrix(c(.6,.4,.8,.2),2,2)
pi_igs_0=abind::abind(pi_gs0,pi_gs0,pi_gs0,along=3)|>
  aperm(c(3,1:2))
epsilon_iba=epsilon_iba_f(3,0.001)


tau_ivgb0=plyr::llply(list(c("aa","cc"),c("ac","ca"),c("gg","tt")),translate_dna_string_vector_to_string_matrix)|>
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


inits=smc_inits(i = i,s=s,v=v,g=g,tau_vgb_0 =tau_vgb_0,pi_igs_0 =pi_igs_0)


desman=smc_custom(n_vsa,
           gs="custom",
           G=2,
           pi_igs_0=pi_igs_0,#set initial value
           block_tau=FALSE,
           bar_epsilon_1=.001,
           alpha_pi=.1,
           n_chains = 3,
           t_max=3,
           smc_kernel=desman_kernel,
           mcmc=TRUE,
           inits=inits)



  
  AM1=smc_custom(n_vsa,
                    gs="custom",
                    G=2,
                    pi_igs_0=pi_igs_0,#set initial value
                    block_tau=TRUE,
                    alpha_tau=NULL,
                    alpha_rho=NULL,
                    alpha_epsilon=NULL,
                    bar_epsilon_1_std=NULL,
                    bar_epsilon_1_mean=NULL,
                    alpha_bar_epsilon=c(1,10),
                    bar_epsilon_1=NULL,
                    kappa_rho=NULL,
                    alpha_pi=.1,
                    n_chains = 2,
                    n_vsa_df=reorder_counts(n_vsa = n_vsa,seed=seed),
                    g=G,
                    t_min=1,
                    t_max,
                    ess_min=NULL,
                    smc_kernel=desman_kernel,
                    trace_all=TRUE,
                    .update_lambda=update_lambda,
                    mcmc=FALSE,
                    inits=NULL,
                    data_tempering=FALSE,
                    bar_epsilon_1_tempering=FALSE,
                    alpha_tau_tempering=FALSE,
                    alpha_tau_seq=NULL,
                    bar_epsilon_seq=NULL,
                    ...){
    
    
    


