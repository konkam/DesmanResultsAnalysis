
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
        einsum::einsum("v,gs,ac->vsgac",1:v,pi_gs,epsilon_ab)+
          einsum::einsum("a,gh,vhb,hs,bc->vsgac",rep(1,4),g_neq_g,tau_vga,pi_gs,epsilon_ba)),
      n_vsa))|>
    plyr::aaply(1:2,function(xi){sample(nucleotides,1,prob = xi)})|>
    translate_dna_matrix_to_binary_array()
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

nu_from_xi<-function(xi){einsum::einsum("vsabg->vsab",xi)}
mu_from_xi<-function(xi){einsum::einsum("vsabg->vsag",xi)}

  
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
  
      rbeta(1,
            c(1,1/delta)+
              c(
                einsum::einsum("vsab,ab->",nu_vsab,a_neq_b),
                einsum::einsum("vsaa->",nu_vsab)))
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
nu_vsab=nu_from_xi(xi=xi_vsabg)
mu_vsab=mu_from_xi(xi=xi_vsabg)
pi_gs=sampler_pi(mu_vsag = mu_vsab,alpha_g = alpha_g)
tau_vga=sampler_tau(n_vsa = n_vsa,tau_vga = tau_vga,pi_gs = pi_gs,epsilon_ba = epsilon_ba)
list(tau_vga=tau_vga,pi_gs=pi_gs,epsilon_ba=epsilon_ba)
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
#'n_vsa_lamba=sub_sample_counts(n_vsa_df,lambda=.5)
#'sub_sample_counts(n_vsa_df,lambda=0)

sub_sample_counts<-function(n_vsa_df,lambda){
  n_vsa_df|>
    dplyr::filter(lambda_threshold<=lambda)|>
    dplyr::select(-lambda_threshold)|>
    plyr::daply(.drop_i=FALSE,
                .drop_o=FALSE,
                .variables=~v+s+a,
                .fun=nrow)}



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
#'tau_vgai=abind::abind(sim1$tau_vga,sim2$tau_vga,along=4)
#'pi_gsi=abind::abind(sim1$pi_gs,sim2$pi_gs,along=3)
#'epsilon_bai=abind::abind(sim1$epsilon_ba,sim2$epsilon_ba,along=3)

smc_logb_prime<-function(n_vsa_lamba,tau_vgai,pi_gsi,epsilon_bai){
  einsum::einsum("vsai,vsa->i",
                 einsum::einsum("vgbi,gsi,bai->vsai",tau_vgai,pi_gsi,epsilon_bai)|>
      log(),
      n_vsa_lambda)}

smc_b_prime<-function(n_vsa_lamba,tau_vgai,pi_gsi,epsilon_bai){
  exp(smc_logb_prime(n_vsa_lamba,tau_vgai,pi_gsi,epsilon_bai))}


#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'plyr::rlply(100,
#'           sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0))->sim
#'n_vsa_df<-reorder_reads(sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)$n_vsa)
#'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
#'tau_vgai=plyr::llply(sim,'[[',"tau_vga")|>c(list(along=4))|>do.call(what=abind::abind)
#'pi_gsi=plyr::llply(sim,'[[',"pi_gs")|>c(list(along=3))|>do.call(what=abind::abind)
#'epsilon_bai=plyr::llply(sim,'[[',"epsilon_ba")|>c(list(along=3))|>do.call(what=abind::abind)
#'smc_ess(n_vsa_lamba,tau_vgai,pi_gsi,epsilon_bai)

smc_ess<-function(n_vsa_lamba,tau_vgai,pi_gsi,epsilon_bai){
  logbs<-smc_logb_prime(n_vsa_lamba,tau_vgai,pi_gsi,epsilon_bai)
  logbs<-logbs-max(logbs)
  sum(exp(logbs))^2/sum(exp(2*logbs))}



#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'plyr::rlply(100,
#'           sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0))->sim
#'n_vsa_df<-reorder_reads(sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)$n_vsa)
#'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
#'tau_vgai=plyr::llply(sim,'[[',"tau_vga")|>c(list(along=4))|>do.call(what=abind::abind)
#'pi_gsi=plyr::llply(sim,'[[',"pi_gs")|>c(list(along=3))|>do.call(what=abind::abind)
#'epsilon_bai=plyr::llply(sim,'[[',"epsilon_ba")|>c(list(along=3))|>do.call(what=abind::abind)
#'ess_min=3
#'max_lambda=1
#'min_lambda.5
#'unaccounted_for=unaccounted_for_f(n_vsa_df,lambda=min_lambda)
#'n_unaccounted_for=nrow(unaccounted_for)
#'smc_new_setup(ess_min=ess_min,
#'     n_vsa_lambda,
#'     unaccounted_for=unaccounted_for,
#'     n_unaccounted_for=n_unaccounted_for,
#'     min_lambda=min_lambda,
#'     max_lambda=1,
#'     tau_vgai,
#'pi_gsi,
#'epsilon_bai)

smc_new_setup<-function(ess_min,
                        n_vsa_lambda,
                        unaccounted_for,
                        n_unaccounted_for=nrow(unaccounted_for),
                        min_lambda,
                        max_lambda=max(unaccounted_for$lambda_threshold),
                        tau_vgai,
                        pi_gsi,
                        epsilon_bai){
  
  print(paste0(min_lambda," _ ",max_lambda))
  if(smc_ess(n_vsa_lamba,tau_vgai,pi_gsi,epsilon_bai)>= ess_min){
    sol=list(n_vsa_lambda=n_vsa_lambda,
         lambda=min_lambda)
    }else{
      if(n_unaccounted_for>0){}
        n=ceiling(n_unaccounted_for/2)
        lambda=max(min_lambda,unaccounted_for[n,]$lambda_threshold)
        n_vsa_lambda_add=sub_sample_counts(unaccounted_for,lambda)
        ess=smc_ess(n_vsa_lamba+n_vsa_lambda_add,tau_vgai,pi_gsi,epsilon_bai)
        if(ess>=ess_min){
          unaccounted_for=unaccounted_for[1:n,]
          n_unaccounted_for=n
          sol=smc_new_setup(ess_min,
                            n_vsa_lambda,
                            unaccounted_for,
                            n_unaccounted_for,
                            min_lambda,
                            max_lambda=lambda,
                            tau_vgai,
                            pi_gsi,
                            epsilon_bai)
      
      
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
                 lambda=max_lambda)
      }else{
      sol=smc_new_setup(ess_min,
                        n_vsa_lambda,
                        unaccounted_for,
                        n_unaccounted_for,
                        min_lambda=lambda,
                        max_lambda=max_lambda,
                        tau_vgai,
                        pi_gsi,
                        epsilon_bai)}
    }}
  
  
  sol
}








#'@description SMC algorithm
#'@param i : integer >0, number of particles
#'@param g : integer >0, number of variants in model
#'@param n_vsa : an array of integers, index by v(position) s (sample) and a (nucleotide a,c,g,t)
#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
#'b=function(lambda){
#'     }
#'  ess_min=1
#'smc_sampler(n_vsa=n_vsa,g=g,t_max=3,i=i,ess_min,b)


smc_sampler<-function(n_vsa,g,t_max,i,ess_min,b){
xi_vsabg=sampler_xi(tau_vga=tau_vga,
                    pi_gs=pi_gs,
                    epsilon_ba = epsilon_ba)
nu_vsab=nu_from_xi(xi=xi_vsabg)
mu_vsab=mu_from_xi(xi=xi_vsabg)
pi_gs=sampler_pi(mu_vsag = mu_vsab,alpha_g = alpha_g)
tau_vga=sampler_tau(n_vsa = n_vsa,tau_vga = tau_vga,pi_gs = pi_gs,epsilon_ba = epsilon_ba)
list(tau_vga=tau_vga,pi_gs=pi_gs,epsilon_ba=epsilon_ba)
}
