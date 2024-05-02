
g_neq_g_f=function(g){1-diag(g)}


#'@description sampler_tau
#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
#'epsilon_ba=.999*diag(4)+.001/3*(1-diag(4))
sampler_tau<-function(tau_vga,pi_gs,epsilon_ba,n_vsa,
                      v=dim(tau_vga)[1],
                      #s=dim(pi_gs)[2],
                      g=dim(pi_gs)[1],
                      g_neq_g=g_neq_g_f(g)){
  exp(
    einsum::einsum(
      equation_string="vsgac,vsc->vga",
      log(
        einsum::einsum("v,gs,ac->vsgac",1:v,pi_gs,epsilon_ba)+
          einsum::einsum("a,gh,vhb,hs,bc->vsgac",rep(1,4),g_neq_g,tau_vga,pi_gs,epsilon_ba)),
      n_vsa))|>
    plyr::aaply(1:2,function(xi){sample(nucleotides,1,prob = xi)})|>
    translate_dna_matrix_to_binary_array()
}


sampler_tau_i<-function(tau_ivga,pi_igs,epsilon_iba,n_vsa,
                        v=dim(tau_ivga)[2],
                        #s=dim(pi_gs)[2],
                        g=dim(pi_igs)[2],
                        g_neq_g=g_neq_g_f(g)){
  exp(
    einsum::einsum(
      equation_string="ivsgac,vsc->ivga",
      log(
        einsum::einsum("v,igs,iac->ivsgac",1:v,pi_igs,epsilon_iab)+
          einsum::einsum("a,gh,ivhb,ihs,ibc->ivsgac",rep(1,4),g_neq_g,tau_ivga,pi_igs,epsilon_iba)),
      n_vsa))|>
    plyr::aaply(1:3,function(xi){sample(nucleotides,1,prob = xi)})|>
    plyr::aaply(1,translate_dna_matrix_to_binary_array)
}

#'@description sampler_mu_nu
#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
#'epsilon_ba=.999*diag(4)+.001/3*(1-diag(4))
#'sampler_xi(tau_vga,pi_gs,epsilon_ba)
sampler_xi<-function(tau_vga,pi_gs,epsilon_ba,
                     v=dim(tau_vga)[1],
                     s=dim(pi_gs)[2],
                     g=dim(pi_gs)[1]){
  ns<-  einsum::einsum("vsa,b,g->vsabg",n_vsa,rep(1,4),rep(1,g))
  ms<-einsum::einsum("vgb,gs,ba->vsabg",tau_vga,pi_gs,epsilon_ba)
  nms=abind::abind(ns,ms,along=6)
  nm=nms[1,1,1,,,]
  nms|>plyr::aaply(1:3,function(nm){
    rmultinom(1,size=nm[1,1,1],prob = nm[,,2])|>
      matrix(4)})
  
}


sampler_xi_i<-function(i,tau_ivga,pi_igs,epsilon_iba,
                       v=dim(tau_vga)[1],
                       s=dim(pi_gs)[2],
                       g=dim(pi_gs)[1]){
  ns<-  einsum::einsum("i,vsa,b,g->vsabg",rep(1,i),n_vsa,rep(1,4),rep(1,g))
  ms<-einsum::einsum("ivgb,igs,iba->ivsabg",tau_ivga,pi_igs,epsilon_iba)
  nms=abind::abind(ns,ms,along=7)
  nm=nms[1,1,1,1,,,]
  nms|>plyr::aaply(1:4,function(nm){
    rmultinom(1,size=nm[1,1,1],prob = nm[,,2])|>
      matrix(4)})
  
}


nu_from_xi<-function(xi){einsum::einsum("vsabg->vsab",xi)}
mu_from_xi<-function(xi){einsum::einsum("vsabg->vsag",xi)}

nu_i_from_xi_i<-function(xi){einsum::einsum("ivsabg->ivsab",xi)}
mu_i_from_xi_i<-function(xi){einsum::einsum("ivsabg->ivsag",xi)}

#'@description sampler_pi
#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
#'xi=sampler_xi(tau_vga,pi_gs,epsilon_ba)
#'mu_vsag=mu_from_xi(xi)
#'alpha_g=rep(.1,g)
#'sampler_pi(mu_vsag,alpha_g)

sampler_pi<-function(mu_vsag,
                     alpha_g){
  
  einsum::einsum("vsag->sg",mu_vsag)|>
    plyr::aaply(1,function(x_g){
      x=dirmult::rdirichlet(1,alpha=alpha_g+x_g)
      x[x<.Machine$double.xmin]<-.Machine$double.xmin
      x/sum(x)})|>
    t()
}


sampler_pi_i<-function(mu_ivsag,
                       alpha_g){
  
  einsum::einsum("ivsag->isg",mu_ivsag)|>
    plyr::aaply(1:2,function(x_g){
      x=dirmult::rdirichlet(1,alpha=alpha_g+x_g)
      x[x<.Machine$double.xmin]<-.Machine$double.xmin
      x/sum(x)})|>
    aperm(c(1,3,2))
}

a_neq_b=g_neq_g_f(4)

#'@description sampler_epsilon_star
#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
#'
#'nu_vsab=nu_from_xi(xi)
#'delta=.1
#'
#'sampler_epsilon_star(nu_vsab,alpha_g)

sampler_epsilon_star<-function(nu_vsab,delta){
  tilde_epsilon=rbeta(1,
        shape1=1+einsum::einsum("vsab,ab->",nu_vsab,a_neq_b),
        shape2=1/delta+einsum::einsum("vsaa->",nu_vsab))
  (1-4/3*tilde_epsilon)*diag(4)+(tilde_epsilon/3)
  }



#'@description SMC kernel
#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
#'





smc_kernel<-function(n_vsa,tau_vga,pi_gs,epsilon_ba,alpha_g,delta){
  xi_vsabg=sampler_xi(tau_vga=tau_vga,
                      pi_gs=pi_gs,
                      epsilon_ba = epsilon_ba)
  nu_vsab<-nu_from_xi(xi=xi_vsabg)
  mu_vsab<-mu_from_xi(xi=xi_vsabg)
  pi_gs<-sampler_pi(mu_vsag = mu_vsab,alpha_g = alpha_g)
  tau_vga<-sampler_tau(n_vsa = n_vsa,tau_vga = tau_vga,pi_gs = pi_gs,epsilon_ba = epsilon_ba)
  epsilon_ba<-sampler_epsilon_star(nu_vsab,delta)
  list(tau_vga=tau_vga,pi_gs=pi_gs,epsilon_ba=epsilon_ba)
}

smc_kernel_i<-function(i,g,n_vsa,tau_ivga,pi_igs,epsilon_iba,alpha_g,delta,g_neq_g =g_neq_g_f(g)){
  xi_ivsabg=sampler_xi_i(i =i, tau_ivga=tau_ivga,
                         pi_igs=pi_igs,
                         epsilon_iba = epsilon_iba)
  #plyr::aaply(nu_ivsab,1,sampler_epsilon_star,delta=delta)
  nu_ivsab=nu_i_from_xi_i(xi=xi_ivsabg)
  mu_ivsab=mu_i_from_xi_i(xi=xi_ivsabg)
  pi_igs=sampler_pi_i(mu_ivsag = mu_ivsab,alpha_g = alpha_g)
  tau_ivga=sampler_tau_i(n_vsa = n_vsa,tau_ivga = tau_ivga,pi_igs = pi_igs,epsilon_iba = epsilon_iba,g_neq_g=g_neq_g )
  list(tau_ivga=tau_ivga,pi_igs=pi_igs,epsilon_iba=epsilon_iba)
}
#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
#'n_vsa_df<-reorder_reads(n_vsa)
#'

reorder_reads<-function(n_vsa,seed=1){
  set.seed(seed)
  n_vsa_df<-n_vsa|>
    plyr::adply(1:3,.fun=function(x){data.frame(V1=rep(1,x))})|>
    dplyr::mutate(rank=rank(runif(dplyr::n())))|>
    dplyr::arrange(rank)|>
    dplyr::mutate(lambda_threshold=dplyr::row_number()/dplyr::n())|>
    dplyr::select(-rank,-V1)
}


#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
#'n_vsa_df<-reorder_reads(n_vsa)
#'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
#'sub_sample_counts(n_vsa_df,lambda=0)

sub_sample_counts<-function(n_vsa_df,lambda){
  n_vsa_df|>
    dplyr::filter(lambda_threshold<=lambda)|>
    dplyr::select(-lambda_threshold)|>
    plyr::daply(.drop_i=FALSE,
                .drop_o=TRUE,
                .variables=~v+s+a,
                .fun=nrow)}
#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1;lambda=.5
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
#'n_vsa_df<-reorder_reads(n_vsa)
#'data_tempering_srs(n_vsa,lambda)

data_tempering_srs<-function(n_vsa,lambda,seed=1,n_vsa_df=reorder_reads(n_vsa,lambda)){
  sub_sample_counts(n_vsa_df,lambda)}


#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1;lambda=.09
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
#'n_vsa_df<-reorder_reads(n_vsa)
#'
#'https://www.census.gov/content/dam/Census/library/working-papers/2014/adrm/rrs2014-07.pdf
reorder_counts<-function(n_vsa,seed=1){
  set.seed(seed)
  H=prod(dim(n_vsa))
  data.frame(n_vsa_v=c(n_vsa),
             i=1:H,
             r=(1:H)[order(runif(H))])
}


#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1;lambda=.09
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
#'n_vsa_df<-reorder_counts(n_vsa,seed=1)
#'data_tempering_stratified(n_vsa,n_vsa_df=n_vsa_df,lambda=lambda)
data_tempering_stratified<-function(n_vsa,n_plus=sum(n_vsa),seed=1,n_vsa_df=reorder_counts(n_vsa,seed=seed),lambda){
  sample_size=ceiling(n_plus*lambda)
  df=n_vsa_df|>
    dplyr::mutate(n_vsa_lambda=n_vsa_v*lambda,
                  allocation0=floor(n_vsa_lambda),
                  left=(n_vsa_v-allocation0/lambda),
                  loss1=-2*lambda+(2*allocation0+1)/n_vsa_v,
                  loss2=1/n_vsa_v)|>
    dplyr::arrange(loss1,loss2,r)
  remaining=sample_size-sum(df$allocation0)
  H=nrow(n_vsa_df)
  while(remaining>0){
    df[1,"allocation0"]<-    df[1,"allocation0"]+1
    loss=df[1,"loss1"]      <-    if( df[1,"allocation0"]>=df[1,"n_vsa_v"]){Inf }else{df[1,"loss2"]}
    i=1
    while(df[i+1,"loss1"]<loss&i<H){
      i=i+1
    }
    #print(i)
    if(i>1){df[1:i,]=df[c(2:i,1),]}
    
  remaining=remaining-1  
    
  }
  
    df$allocation0[df$i]|>array(dim = dim(n_vsa),dimnames=dimnames(n_vsa))
}
  
  
  #'@examples
  #'n=1000;v=20;g=5;s=3;alpha0=.1
  #'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
  #'n_vsa_df<-reorder_reads(n_vsa)
  #'unaccounted_for=unaccounted_for_f(n_vsa_df,lambda=.5)
  
  unaccounted_for_f<-function(n_vsa_df,lambda){
    n_vsa_df|>dplyr::filter(lambda_threshold>lambda)}
  
  
  #'@examples
  #'n=1000;v=20;g=5;s=3;alpha0=.1
  #'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)->sim1
  #'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)->sim2
  #'n_vsa_df<-reorder_reads(sim1$n_vsa)
  #'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
  #'tau_ivga=abind::abind(sim1$tau_vga,sim2$tau_vga,along=4)|>aperm(c(4,1:3))
  #'pi_igs=abind::abind(sim1$pi_gs,sim2$pi_gs,along=3)|>aperm(c(3,1:2))
  #'epsilon_iba=abind::abind(sim1$epsilon_ba,sim2$epsilon_ba,along=3)|>aperm(c(3,1:2))
  #'smc_logb_prime(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
  
  smc_logb_prime<-function(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba){
    einsum::einsum("ivsa,vsa->i",
                   einsum::einsum("ivgb,igs,iba->ivsa",tau_ivga,pi_igs,epsilon_iba)|>
                     log(),
                   n_vsa_lambda)}
  
  smc_b_prime<-function(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba){
    exp(smc_logb_prime(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba))}
  
  
  #'@examples
  #'n=1000;v=20;g=5;s=3;alpha0=.1
  #'plyr::rlply(100,
  #'           sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0))->sim
  #'n_vsa_df<-reorder_reads(sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)$n_vsa)
  #'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
  #'tau_ivga=plyr::llply(sim,'[[',"tau_vga")|>c(list(along=4))|>do.call(what=abind::abind)
  #'pi_igs=plyr::llply(sim,'[[',"pi_gs")|>c(list(along=3))|>do.call(what=abind::abind)
  #'epsilon_iba=plyr::llply(sim,'[[',"epsilon_ba")|>c(list(along=3))|>do.call(what=abind::abind)
  #'smc_ess(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
  
  smc_ess<-function(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba){
    logbs<-smc_logb_prime(n_vsa_lambda=n_vsa_lambda,tau_ivga=tau_ivga,pi_igs=pi_igs,epsilon_iba=epsilon_iba)
    logbs<-logbs-max(logbs)
    sum(exp(logbs))^2/sum(exp(2*logbs))}
  
  
  
  #'@examples
  #'n=1000;v=20;g=5;s=3;alpha0=.1
  #'plyr::rlply(100,
  #'           sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0))->sim
  #'           n_vsa=sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)$n_vsa
  #'n_vsa_df<-reorder_reads(n_vsa)
  #'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
  #'tau_ivga=plyr::llply(sim,'[[',"tau_vga")|>c(list(along=4))|>do.call(what=abind::abind)
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
  #'     tau_ivga,
  #'pi_igs,
  #'epsilon_iba)
  
  smc_new_setup<-function(ess_min,
                          n_vsa_lambda,
                          unaccounted_for,
                          n_unaccounted_for=nrow(unaccounted_for),
                          min_lambda,
                          max_lambda=max(unaccounted_for$lambda_threshold),
                          tau_ivga,
                          pi_igs,
                          epsilon_iba){
    
    #print(paste0(min_lambda," _ ",max_lambda))
    ess=smc_ess(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
    if(ess>= ess_min|n_unaccounted_for==0){
      sol=list(n_vsa_lambda=n_vsa_lambda,
               lambda=min_lambda,
               ess=ess)
    }else{
      n=ceiling(n_unaccounted_for/2)
      lambda=max(min_lambda,unaccounted_for[n,]$lambda_threshold)
      n_vsa_lambda_add=sub_sample_counts(unaccounted_for,lambda)
      ess=smc_ess(n_vsa_lambda+n_vsa_lambda_add,tau_ivga,pi_igs,epsilon_iba)
      if(ess>=ess_min){
        unaccounted_for=unaccounted_for[1:n,]
        n_unaccounted_for=n
        sol=smc_new_setup(ess_min,
                          n_vsa_lambda,
                          unaccounted_for,
                          n_unaccounted_for,
                          min_lambda,
                          max_lambda=lambda,
                          tau_ivga,
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
                            tau_ivga,
                            pi_igs,
                            epsilon_iba)}
      }}
    
    
    sol
  }
  
  
  #'@examples
  #'n=1000;v=20;g=5;s=3;alpha0=.1
  #'plyr::rlply(100,
  #'           sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0))->sim
  #'n_vsa=sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)$n_vsa
  #'tau_ivga=plyr::llply(sim,'[[',"tau_vga")|>c(list(along=4))|>do.call(what=abind::abind)
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
  #'tau_ivga=tau_ivga,
  #'pi_igs=pi_igs,
  #'  epsilon_iba=epsilon_iba)
  
  update_lambda<-function(ess_min,
                          n_vsa,
                          n_plus=sum(n_vsa),
                          seed=1,
                          n_vsa_df=reorder_counts(n_vsa,seed=seed),
                          min_samplesize=0,
                          max_samplesize=n_plus,
                          tau_ivga,
                          pi_igs,
                          epsilon_iba){
    
    if(min_samplesize==max_samplesize){
      lambda=min_samplesize/n_plus
      n_vsa_lambda=data_tempering_stratified(n_vsa=n_vsa,n_plus=n_plus,n_vsa_df=n_vsa_df,lambda=lambda)
      ess=smc_ess(n_vsa_lambda=n_vsa_lambda,tau_ivga=tau_ivga,pi_igs=pi_igs,epsilon_iba=epsilon_iba)
      sol=list(n_vsa_lambda=n_vsa_lambda,
               lambda=lambda,ess=ess)
    }else{
      samplesize=floor((min_samplesize+max_samplesize)/2)
      lambda=samplesize/n_plus
      n_vsa_lambda=data_tempering_stratified(n_vsa=n_vsa,n_plus=n_plus,n_vsa_df=n_vsa_df,lambda=lambda)
      ess=smc_ess(n_vsa_lambda=n_vsa_lambda,tau_ivga=tau_ivga,pi_igs=pi_igs,epsilon_iba=epsilon_iba)
      
      #print(paste0(min_samplesize," _ ",samplesize," _ ",max_samplesize," _ ",ess))
      
      if(ess>=ess_min){
        sol<-update_lambda(ess_min=ess_min,
                       n_vsa=n_vsa,
                       n_plus=n_plus,
                       n_vsa_df=n_vsa_df,
                       min_samplesize=min_samplesize,
                       max_samplesize=samplesize,
                       tau_ivga=tau_ivga,
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
                           tau_ivga=tau_ivga,
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
  #'n=1000;v=20;g=5;s=3;alpha0=.1
  #'n_vsa=sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)$n_vsa
  #'n_vsa_df<-reorder_reads(n_vsa)
  #'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
  #'i=10
  #'delta=1e-3
  #'ess_min=3
  #'max_lambda=1
  #'min_lambda=.5
  #'unaccounted_for=unaccounted_for_f(n_vsa_df,lambda=min_lambda)
  #' t_max=30
  #' t_min=30
  #' n_plus=sum(n_vsa)
  #' new=smc_sampler(n_vsa=n_vsa,g=g,t_max=t_max,i=i,ess_min=ess_min)
  
  
  smc_sampler<-function(n_vsa,
                        seed=1,
                        n_plus=sum(n_vsa),
                        n_vsa_df=reorder_counts(n_vsa = n_vsa,seed=seed),
                        g,
                        i,
                        t_min=1,
                        t_max,
                        ess_min,
                        delta=1e-3,
                        alpha0=.1,
                        trace_all=TRUE,
                        .update_lambda=update_lambda){
    
    lambda=0
    alpha_g=rep(alpha0,g)
    t=0
    v=dim(n_vsa)[1]
    s=dim(n_vsa)[2]
    epsilon_iba=rbeta(i,1,1/delta)|>plyr::aaply(1,epsilon_ba_f)
    tau_ivga <- plyr::raply(i,sim_tau_vga(v = v, g = g))
    pi_igs <- plyr::raply(i,sim_pi_gs(g = g, s = s, alpha0 = alpha0))
    #w<-smc_b_prime(0*n_vsa,tau_ivga,pi_igs,epsilon_iba)
    new=list(epsilon_iba=epsilon_iba,tau_ivga=tau_ivga,pi_igs=pi_igs)
    ww<-rep(1/i,i)#w/sum(w)
    if(trace_all){
      traces=list(epsilon_iba=array(epsilon_iba,c(1,dim(epsilon_iba))),
                 tau_ivga=array(tau_ivga,c(1,dim(tau_ivga))),
                 pi_igs=array(pi_igs,c(1,dim(pi_igs))),
                 ww=array(ww,c(1,dim(ww))),
                 sample_i=array(1:i,c(1,i)),
                 ess=1)
                 
    }
    while(((lambda<1)|(t<=min(t_max,t_min)))&(t<max(t_max,t_min))){
      print(paste0("t : ",t," || lambda :",lambda))
      t=t+1
      n_vsa_lambda=data_tempering_stratified(n_vsa=n_vsa,n_plus=n_plus,n_vsa_df=n_vsa_df,lambda=lambda)
      sample_i=sample(x=i,size=i,replace=(i>1),prob=ww)
      new=plyr::llply(new,resample_array,dimension=1,selection=sample_i)
      
      new=with(new,
               plyr::alply(seq_len(i),1,
                      function(i){
                        smc_kernel(n_vsa_lambda,
                                   tau_ivga[i,,,],
                                   pi_igs[i,,],
                                   epsilon_iba[i,,],
                                   alpha_g,delta)}))|>
        (function(x){plyr::alply(c("tau_vga","pi_gs","epsilon_ba"),1,function(parameter){
          plyr::laply(x,`[[`,parameter)})})()|>
        setNames(c("tau_ivga","pi_igs","epsilon_iba"))
      
          
      w<-do.call(what=smc_b_prime,c(list(n_vsa_lambda),new[c("tau_ivga","pi_igs","epsilon_iba")]))
      ww<-w/sum(w)
      lambda_n_vsa=.update_lambda(ess_min,
                    n_vsa=n_vsa,
                    n_plus=n_plus,seed=seed,n_vsa_df=n_vsa_df,min_samplesize=0,
                    max_samplesize=n_plus,
                    tau_ivga,
                    pi_igs,
                    epsilon_iba)
      lambda=lambda_n_vsa$lambda
      n_vsa_lambda=lambda_n_vsa$n_vsa_lambda
      
      
      
      if(trace_all){
        traces=list(epsilon_iba=abind::abind(traces$epsilon_iba,new$epsilon_iba,along=1),
                    pi_igs=abind::abind(traces$pi_igs,new$pi_igs,along=1),
                    epsilon_iba=abind::abind(traces$epsilon_iba,new$epsilon_iba,along=1),
                    ww=abind::abind(traces$ww,array(ww,c(1,i)),along=1),
                    epsilon_iba=abind::abind(traces$epsilon_iba,new$epsilon_iba,along=1),
                    sample_i=abind::abind(traces$sample_i,array(sample_i,c(1,i)),along=1),
                    ess=c(traces$ess,(sum(ww))^2/(sum(ww^2))))}
        
      
    }
    list(sample=new,traces=if(trace_all){traces}else{NULL})
    
  }

  