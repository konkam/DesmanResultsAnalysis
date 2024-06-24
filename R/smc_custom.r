

#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)->sim1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)->sim2
#'n_vsa_df<-reorder_reads(sim1$n_vsa)
#'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
#'tau_ivgb=abind::abind(sim1$tau_vgb,sim2$tau_vgb,along=4)|>aperm(c(4,1:3))
#'pi_igs=abind::abind(sim1$pi_gs,sim2$pi_gs,along=3)|>aperm(c(3,1:2))
#'epsilon_iba=abind::abind(sim1$epsilon_ba,sim2$epsilon_ba,along=3)|>aperm(c(3,1:2))
#'smc_logb_prime(n_vsa_lambda,tau_ivgb,pi_igs,epsilon_iba)

smc_logb_prime<-function(n_vsa,tau_ivgb,pi_igs,epsilon_iba){
  einsum::einsum("ivsa,vsa->i",
                 einsum::einsum("ivgb,igs,iba->ivsa",tau_ivgb,pi_igs,epsilon_iba)|>
                   log(),
                 n_vsa)}

smc_b_prime<-function(n_vsa,tau_ivgb,pi_igs,epsilon_iba){
  exp(smc_logb_prime(n_vsa,tau_ivgb,pi_igs,epsilon_iba))}


#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'plyr::rlply(100,
#'           sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi))->sim
#'n_vsa_df<-reorder_reads(sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)$n_vsa)
#'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
#'tau_ivgb=plyr::llply(sim,'[[',"tau_vgb")|>c(list(along=4))|>do.call(what=abind::abind)
#'pi_igs=plyr::llply(sim,'[[',"pi_gs")|>c(list(along=3))|>do.call(what=abind::abind)
#'epsilon_iba=plyr::llply(sim,'[[',"epsilon_ba")|>c(list(along=3))|>do.call(what=abind::abind)
#'smc_ess(n_vsa_lambda,tau_ivgb,pi_igs,epsilon_iba)

smc_ess<-function(n_vsa_lambda,tau_ivgb,pi_igs,epsilon_iba){
  logbs<-smc_logb_prime(n_vsa_lambda=n_vsa_lambda,tau_ivgb=tau_ivgb,pi_igs=pi_igs,epsilon_iba=epsilon_iba)
  logbs<-logbs-max(logbs)
  sum(exp(logbs))^2/sum(exp(2*logbs))}



#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'plyr::rlply(100,
#'           sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi))->sim
#'           n_vsa=sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)$n_vsa
#'n_vsa_df<-reorder_reads(n_vsa)
#'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
#'tau_ivgb=plyr::llply(sim,'[[',"tau_vgb")|>c(list(along=4))|>do.call(what=abind::abind)
#'pi_igs=plyr::llply(sim,'[[',"pi_gs")|>c(list(along=3))|>do.call(what=abind::abind)
#'epsilon_iba=plyr::llply(sim,'[[',"epsilon_ba")|>c(list(along=3))|>do.call(what=abind::abind)
#'ess_min=3
#'max_lambda=1
#'min_lambda=.5
#'unaccounted_for=unaccounted_for_f(n_vsa_df,lambda=min_lambda)
#'n_unaccounted_for=nrow(unaccounted_for)
#'smc_new_setup(ess_min=ess_min,
#'     n_vsa_lambda,
#'     unaccounted_for=unaccounted_for,
#'     n_unaccounted_for=n_unaccounted_for,
#'     min_lambda=min_lambda,
#'     max_lambda=1,
#'     tau_ivgb,
#'pi_igs,
#'epsilon_iba)

smc_new_setup<-function(ess_min,
                        n_vsa_lambda,
                        unaccounted_for,
                        n_unaccounted_for=nrow(unaccounted_for),
                        min_lambda,
                        max_lambda=max(unaccounted_for$lambda_threshold),
                        tau_ivgb,
                        pi_igs,
                        epsilon_iba){
  
  #print(paste0(min_lambda," _ ",max_lambda))
  ess=smc_ess(n_vsa_lambda,tau_ivgb,pi_igs,epsilon_iba)
  if(ess>= ess_min|n_unaccounted_for==0){
    sol=list(n_vsa_lambda=n_vsa_lambda,
             lambda=min_lambda,
             ess=ess)
  }else{
    n=ceiling(n_unaccounted_for/2)
    lambda=max(min_lambda,unaccounted_for[n,]$lambda_threshold)
    n_vsa_lambda_add=sub_sample_counts(unaccounted_for,lambda)
    ess=smc_ess(n_vsa_lambda+n_vsa_lambda_add,tau_ivgb,pi_igs,epsilon_iba)
    if(ess>=ess_min){
      unaccounted_for=unaccounted_for[1:n,]
      n_unaccounted_for=n
      sol=smc_new_setup(ess_min,
                        n_vsa_lambda,
                        unaccounted_for,
                        n_unaccounted_for,
                        min_lambda,
                        max_lambda=lambda,
                        tau_ivgb,
                        pi_igs,
                        epsilon_iba)
      
      
    }else{
      n_vsa_lambda<-n_vsa_lambda+
        unaccounted_for[1:n,]|>
        dplyr::select(-lambda_threshold)|>
        plyr::daply(.drop_i=FALSE,
                    .drop_o=TRUE,
                    .variables=~v+s+a,
                    .fun=nrow)
      unaccounted_for=unaccounted_for[-(1:n),]
      n_unaccounted_for=n_unaccounted_for-n
      
      if(n==0){
        sol=list(n_vsa_lambda=n_vsa_lambda,
                 lambda=max_lambda,ess=ess)
      }else{
        sol=smc_new_setup(ess_min,
                          n_vsa_lambda,
                          unaccounted_for,
                          n_unaccounted_for,
                          min_lambda=lambda,
                          max_lambda=max_lambda,
                          tau_ivgb,
                          pi_igs,
                          epsilon_iba)}
    }}
  
  
  sol
}


#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'plyr::rlply(100,
#'           sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi))->sim
#'n_vsa=sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)$n_vsa
#'tau_ivgb=plyr::llply(sim,'[[',"tau_vgb")|>c(list(along=4))|>do.call(what=abind::abind)
#'pi_igs=plyr::llply(sim,'[[',"pi_gs")|>c(list(along=3))|>do.call(what=abind::abind)
#'epsilon_iba=plyr::llply(sim,'[[',"epsilon_ba")|>c(list(along=3))|>do.call(what=abind::abind)
#'ess_min=3
#'max_lambda=1
#'min_samplesize=.5*n_plus
#'max_samplesize=n_plus
#'n_plus=sum(n_vsa)
#'n_vsa_df=reorder_counts(n_vsa,seed=1)
#'update_lambda(ess_min=ess_min,
#'n_vsa=n_vsa,
#'n_vsa_df=n_vsa_df,
#'n_plus=n_plus,
#'min_samplesize=min_samplesize,
#'max_samplesize=max_samplesize,
#'tau_ivgb=tau_ivgb,
#'pi_igs=pi_igs,
#'  epsilon_iba=epsilon_iba)

update_lambda<-function(ess_min,
                        n_vsa,
                        n_plus=sum(n_vsa),
                        seed=1,
                        n_vsa_df=reorder_counts(n_vsa,seed=seed),
                        min_samplesize=0,
                        max_samplesize=n_plus,
                        tau_ivgb,
                        pi_igs,
                        epsilon_iba){
  
  if(min_samplesize==max_samplesize){
    lambda=min_samplesize/n_plus
    n_vsa_lambda=data_tempering_stratified(n_vsa=n_vsa,n_plus=n_plus,n_vsa_df=n_vsa_df,lambda=lambda)
    ess=smc_ess(n_vsa_lambda=n_vsa_lambda,tau_ivgb=tau_ivgb,pi_igs=pi_igs,epsilon_iba=epsilon_iba)
    sol=list(n_vsa_lambda=n_vsa_lambda,
             lambda=lambda,ess=ess)
  }else{
    samplesize=floor((min_samplesize+max_samplesize)/2)
    lambda=samplesize/n_plus
    n_vsa_lambda=data_tempering_stratified(n_vsa=n_vsa,n_plus=n_plus,n_vsa_df=n_vsa_df,lambda=lambda)
    ess=smc_ess(n_vsa_lambda=n_vsa_lambda,tau_ivgb=tau_ivgb,pi_igs=pi_igs,epsilon_iba=epsilon_iba)
    
    #print(paste0(min_samplesize," _ ",samplesize," _ ",max_samplesize," _ ",ess))
    
    if(ess>=ess_min){
      sol<-update_lambda(ess_min=ess_min,
                         n_vsa=n_vsa,
                         n_plus=n_plus,
                         n_vsa_df=n_vsa_df,
                         min_samplesize=min_samplesize,
                         max_samplesize=samplesize,
                         tau_ivgb=tau_ivgb,
                         pi_igs=pi_igs,
                         epsilon_iba=epsilon_iba)
    }else{
      if(max_samplesize-samplesize==1){
        samplesize<-max_samplesize
      }
      sol<-update_lambda(ess_min=ess_min,
                         n_vsa=n_vsa,
                         n_plus=n_plus,
                         n_vsa_df=n_vsa_df,
                         min_samplesize=samplesize,
                         max_samplesize=max_samplesize,
                         tau_ivgb=tau_ivgb,
                         pi_igs=pi_igs,
                         epsilon_iba=epsilon_iba)}}
  
  
  sol
}

#'@examples 
#'the_array=array(rnorm(36),c(3,4,3))
#'dimension=2
#'selection=3:4
#'resample_array(the_array,dimension,selection)|>dim()
#'identical(resample_array(the_array,dimension,selection),the_array[,3:4,])
resample_array<-function(the_array,dimension=1,selection=1:(dim(the_array)[dimension])){
  a=plyr::alply(dim(the_array),1,seq_len)
  a[[dimension]]<-selection#a[dimension]<-if(is.vector(selection)){list(selection)}else{selection}
  do.call(what=`[`,c(list(the_array),a,list(drop=FALSE)))
}




#'@description SMC algorithm
#'@param i : integer >0, number of particles
#'@param g : integer >0, number of variants in model
#'@param n_vsa : an array of integers, index by v(position) s (sample) and a (nucleotide a,c,g,t)
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim=sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)
#'n_vsa=sim$n_vsa
#'i=100
#'init=list(pi_igs=plyr::raply(i,sim$pi_gs),
#'epsilon_iba=plyr::raply(i,sim$epsilon_ba),
#'tau_ivgb=plyr::raply(i,sim$tau_vgb))
#'n_vsa_df=reorder_counts(n_vsa = n_vsa,seed=seed)
#'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
#'delta=1e-3
#'ess_min=60
#'max_lambda=1
#'min_lambda=.5
#'unaccounted_for=unaccounted_for_f(n_vsa_df,lambda=min_lambda)
#' t_max=30
#' t_min=30
#' n_plus=sum(n_vsa)
#' trace_all=TRUE
#' .update_lambda=update_lambda
#' init=NULL 
#' shape_epsilon=c(1,1000)
#' mcmc=TRUE
#' new=smc_sampler(n_vsa=n_vsa,g=g,t_min=t_min,t_max=t_max,i=i,ess_min=ess_min,init=NULL)
#' save(new,file="new.rda")
#' new$traces$ess
#' 
smc_sampler<-function(n_vsa,
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
                      mcmc=FALSE){
  
  lambda=0.1
  rep_alpha_pi=rep(alpha_pi,g)
  t=0
  v=dim(n_vsa)[1]
  s=dim(n_vsa)[2]
  new=if(is.null(init)){{
    epsilon_iba=rbeta(i,shape_epsilon[1],shape_epsilon[2])|>plyr::aaply(1,epsilon_ba_f)
    tau_ivgb <- plyr::raply(i,sim_tau_vgb(v = v, g = g))
    pi_igs <- plyr::raply(i,sim_pi_gs(g = g, s = s, alpha_pi = alpha_pi))
    list(epsilon_iba=epsilon_iba,tau_ivgb=tau_ivgb,pi_igs=pi_igs)}}else{init}
  
  #w<-smc_b_prime(0*n_vsa,tau_ivgb,pi_igs,epsilon_iba)
  if(!mcmc){ww<-rep(1/i,i)}#w/sum(w)}
  if(trace_all){
    traces=list(epsilon_iba=array(new$epsilon_iba,c(1,dim(new$epsilon_iba))),
                tau_ivgb=array(new$tau_ivgb,c(1,dim(new$tau_ivgb))),
                pi_igs=array(new$pi_igs,c(1,dim(new$pi_igs))),
                ww=array(ww,c(1,dim(as.array(ww)))),
                sample_i=array(1:i,c(1,i)),
                ess=1)
    
  }
  if(mcmc){n_vsa_lambda=n_vsa;lambda=0;sample_i=1:i}
  while(((lambda<1)|(t<=min(t_max,t_min)))&(t<max(t_max,t_min))){
    print(paste0("t : ",t," || lambda :",lambda))
    t=t+1
    if(!mcmc){
      n_vsa_lambda=data_tempering_stratified(n_vsa=n_vsa,n_plus=n_plus,n_vsa_df=n_vsa_df,lambda=lambda)
      sample_i=sample(x=i,size=i,replace=(i>1),prob=ww)}
    
    new0=new
    new=plyr::llply(new,resample_array,dimension=1,selection=sample_i)
    new=with(new,
             plyr::alply(seq_len(i),1,
                         function(i){
                           smc_kernel(n_vsa_lambda,
                                      tau_ivgb[i,,,],
                                      pi_igs[i,,],
                                      epsilon_iba[i,,],
                                      rep_alpha_pi,
                                      delta)}))|>
      (function(x){plyr::alply(c("tau_vgb","pi_gs","epsilon_ba"),1,function(parameter){
        plyr::laply(x,`[[`,parameter)})})()|>
      setNames(c("tau_ivgb","pi_igs","epsilon_iba"))
    
    if(!mcmc){
      w<-do.call(what=smc_b_prime,c(list(n_vsa_lambda=n_vsa_lambda),new[c("tau_ivgb","pi_igs","epsilon_iba")]))
      w<-w/max(w)
      ww<-w/sum(w)
      ess=1/sum(ww^2)
      sample_size=sum(n_vsa_lambda)
      lambda_n_vsa=.update_lambda(ess_min,
                                  n_vsa=n_vsa,
                                  n_plus=n_plus,seed=seed,
                                  n_vsa_df=n_vsa_df,
                                  min_samplesize=sample_size,
                                  max_samplesize=n_plus,
                                  tau_ivgb,
                                  pi_igs,
                                  epsilon_iba)
      lambda=max(lambda_n_vsa$lambda,lambda*1.02)
      n_vsa_lambda=lambda_n_vsa$n_vsa_lambda}
    
    
    
    if(trace_all){
      traces=list(epsilon_iba=abind::abind(traces$epsilon_iba,new$epsilon_iba,along=1),
                  pi_igs=abind::abind(traces$pi_igs,new$pi_igs,along=1),
                  tau_ivgb=abind::abind(traces$tau_ivgb,new$tau_ivgb,along=1),
                  ww=abind::abind(traces$ww,array(ww,c(1,i)),along=1),
                  sample_i=abind::abind(traces$sample_i,array(sample_i,c(1,i)),along=1),
                  ess=c(traces$ess,1/(sum(ww^2))))}
    
    
  }
  list(sample=new,traces=if(trace_all){traces}else{NULL})
}








#'@description SMC algorithm
#'@param i : integer >0, number of particles
#'@param g : integer >0, number of variants in model
#'@param n_vsa : an array of integers, index by v(position) s (sample) and a (nucleotide a,c,g,t)
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1;i=7
#'sim=sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)
#'n_vsa=sim$n_vsa
#'i=3
#'  inits=smc_inits(v=20,
#'             s=s,
#'             g=g,
#'             i=i,
#'             tau_vgb_0=NULL,
#'             pi_gs_0=NULL,
#'             gs="jags",fixed_tau=FALSE,alpha_tau=NULL,alpha_epsilon=NULL,
#'             alpha_bar_epsilon=c(1,100),kappa_rho=NULL,alpha_rho=NULL,
#'             bar_epsilon_1_std=NULL,bar_epsilon_1_mean=NULL,alpha_pi=.1) 
#'seed=1
#'n_vsa_df=reorder_counts(n_vsa = n_vsa,seed=seed)
#'#n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
#'delta=1e-3
#'ess_min=60
#'max_lambda=1
#'min_lambda=.5
#'#unaccounted_for=unaccounted_for_f(n_vsa_df,lambda=min_lambda)
#' t_max=30
#' t_min=30
#' n_plus=sum(n_vsa)
#' trace_all=TRUE
#' .update_lambda=update_lambda
#' shape_epsilon=c(1,1000)
#' mcmc=TRUE
#'n_chains=300 
#'alpha_tau=.1
#'kappa_rho=NULL
#'alpha_rho=NULL
#'pi_igs_0=NULL
#'pi_gs_0=NULL
#'alpha_epsilon=NULL
#'block_tau=FALSE
#'tau_vgb_0=NULL
#'alpha_epsilon=NULL
#'bar_epsilon_1_std=NULL
#'bar_epsilon_1_mean=NULL
#'alpha_bar_epsilon=NULL#c(1,10)
#'bar_epsilon_1=.1
#'data_tempering=FALSE
#'bar_epsilon_1_tempering=FALSE
#'alpha_tau_tempering=FALSE
#'#'tau_vgb=NULL
#'alpha_tau_seq=NULL
#'bar_epsilon_1_seq=(function(x){.001+0.249*(1/x)^{1/4}})(1:1000)



smc_custom<-function(n_vsa,
                     gs="jags",
                     tau_vgb=NULL,
                     G=if(!is.null(tau_vgb)){dim(tau_vgb)[2]}else{5},
                     pi_igs_0=NULL,#set initial value
                     pi_gs_0=NULL,#set initial value
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

  fixed_bar_epsilon=!is.null(bar_epsilon_1)
  fixed_tau=!is.null(tau_vgb)
  relax_tau=!is.null(alpha_tau)
  relax_rho=!is.null(kappa_rho)|!is.null(alpha_rho)
  fixed_alpha_rho=!is.null(alpha_rho)
  constrained_epsilon_matrix=is.null(alpha_epsilon)
  if(is.null(inits)){
    inits=smc_inits(v=v,
                    s=s,
                    g=g,
                    i=n_chains,
                    gs=gs,
                    tau_vgb_0=tau_vgb_0,#set initial value
                    tau_ivgb_0=tau_vgb_0,#set initial value
                    fixed_tau=fixed_tau,
                    alpha_tau=alpha_tau,
                    pi_igs_0=pi_igs_0,#set initial value
                    pi_gs_0=pi_gs_0,#set initial value
                    fixed_pi=fixed_pi,
                    alpha_pi=alpha_pi,
                  bar_epsilon_1_0=bar_epsilon_1_0,
                  alpha_bar_epsilon=alpha_bar_epsilon,
                  bar_epsilon_1_std=bar_epsilon_1_std,
                  bar_epsilon_1_mean=bar_epsilon_1_mean,
                  fixed_bar_epsilon=fixed_bar_epsilon,
                  epsilon_ba_0=epsilon_ba_0,
                  alpha_epsilon=alpha_epsilon,
                  kappa_rho=kappa_rho,
                  alpha_rho=alpha_rho)
    }
  
  smc_kernel=kernel_f(fixed_bar_epsilon=fixed_bar_epsilon,
                      constrained_epsilon_matrix=constrained_epsilon_matrix,
                      block_tau=block_tau,
                      fixed_tau=fixed_tau,
                      relax_tau=relax_tau,
                      relax_rho=relax_rho)
  
  
  param=smc_fixed_f(n_vsa=n_vsa,
             gs="custom",
             G=g,
             alpha_pi=alpha_pi,
             alpha_epsilon=alpha_epsilon,
             bar_epsilon_1_std=bar_epsilon_1_std,
             bar_epsilon_1_mean=bar_epsilon_1_mean,
             alpha_bar_epsilon=alpha_bar_epsilon,
             bar_epsilon_1=bar_epsilon_1,
             tau_vgb=tau_vgb,
             alpha_tau=alpha_tau,
             kappa_rho=kappa_rho,
             alpha_rho=alpha_rho,
             i=n_chains)
  
  theta=inits
  fixed=param
  
  
    
  lambda=0.1
  t=0
  n_plus=sum(n_vsa)
  v=dim(n_vsa)[1]
  s=dim(n_vsa)[2]
 
  
  #w<-smc_b_prime(0*n_vsa,tau_ivgb,pi_igs,epsilon_iba)
  if(!mcmc){ww<-rep(1/n_chains,n_chains)}#w/sum(w)}
  if(trace_all){
    traces=list()
    if(is.element("epsilon_iba",names(theta))){
      traces=c(traces,list(epsilon_tiba=array(theta$epsilon_iba,c(1,dim(theta$epsilon_iba)))))}
    if(is.element("pi_igs",names(theta))){
      traces=c(traces,list(pi_tigs=array(theta$pi_igs,c(1,dim(theta$pi_igs)))))}
    if(is.element("tau_ivgb",names(theta))){
      traces=c(traces,list(tau_tivgb=array(theta$tau_ivgb,c(1,dim(theta$tau_ivgb)))))}
    traces=c(traces,list(
            ww_ti=array(ww,c(1,dim(as.array(ww)))),
                sample_ti=array(1:n_chains,c(1,n_chains)),
                ess_t=1))
    
  }
  sample_i=1:n_chains
  if(!data_tempering){n_vsa_lambda=n_vsa;lambda=0}
  while(((lambda<1)|(t<=min(t_max,t_min)))&(t<max(t_max,t_min))){
    print(paste0("t : ",t," || lambda :",lambda))
    t=t+1
    if(data_tempering){
      n_vsa_lambda=data_tempering_stratified(n_vsa=n_vsa,n_plus=n_plus,n_vsa_df=n_vsa_df,lambda=lambda)
      fixed$n_vsa=n_vsa_lambda
      }
    
    if(data_tempering){
      n_vsa_lambda=data_tempering_stratified(n_vsa=n_vsa,n_plus=n_plus,n_vsa_df=n_vsa_df,lambda=lambda)
      fixed$n_vsa=n_vsa_lambda
    }
    
    if(bar_epsilon_1_tempering){
      fixed$bar_epsilon_1=bar_epsilon_1_seq[min(t,length(bar_epsilon_1_seq))]
      fixed$bar_epsilon=c(fixed$bar_epsilon_1,1-fixed$bar_epsilon_1)
      fixed$epsilon_ba=epsilon_ba_f(bar_epsilon_1=fixed$bar_epsilon_1)
      epsilon_iba=epsilon_iba_f(i=n_chains,bar_epsilon_1=fixed$bar_epsilon_1)
      
      
    }
      if(!mcmc&t>1){
      sample_i=sample(x=n_chains,size=n_chains,replace=(n_chains>1),prob=ww)
      }
    
    obs=list(n_vsa=n_vsa)
    if(!mcmc){theta=plyr::llply(theta,resample_array,dimension=1,selection=sample_i)}
    theta=smc_kernel(theta = theta,
                     fixed=fixed)
    
    if(!mcmc){
      w<-do.call(what=smc_logb_prime,c(theta,fixed)[c("n_vsa","tau_ivgb","pi_igs","epsilon_iba")])
      w<-exp(w-max(w))
      ww<-w/sum(w)
      ess=1/sum(ww^2)
      if(data_tempering){
      sample_size=sum(n_vsa_lambda)
      lambda_n_vsa=.update_lambda(ess_min,
                                  n_vsa=n_vsa,
                                  n_plus=n_plus,seed=seed,
                                  n_vsa_df=n_vsa_df,
                                  min_samplesize=sample_size,
                                  max_samplesize=n_plus,
                                  tau_ivgb,
                                  pi_igs,
                                  epsilon_iba)
      lambda=max(lambda_n_vsa$lambda,lambda*1.02)
      n_vsa_lambda=lambda_n_vsa$n_vsa_lambda}
      
      }
    
    
    
    if(trace_all){
      if(is.element("epsilon_iba",names(theta))){
        traces$epsilon_tiba=abind::abind(traces$epsilon_tiba,theta$epsilon_iba,along=1)}
      if(is.element("pi_igs",names(theta))){
        traces$pi_tigs=abind::abind(traces$pi_tigs,theta$pi_gs,along=1)}
      if(is.element("tau_ivgb",names(theta))){
        traces$tau_tivgb=abind::abind(traces$tau_tivgb,theta$tau_ivgb,along=1)}
      traces$ww_ti=abind::abind(traces$ww_ti,array(ww,c(1,n_chains)),along=1)
      traces$sample_ti=abind::abind(traces$sample_ti,array(sample_i,c(1,n_chains)),along=1)
      traces$ess_t=c(traces$ess_t,1/(sum(ww^2)))
      }
  }
  list(theta=theta,traces=if(trace_all){traces}else{NULL})
}

