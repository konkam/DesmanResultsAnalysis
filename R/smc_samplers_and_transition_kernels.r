g_neq_g_f=function(g){1-diag(g)}


#'@description sampler_tau
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha_pi=alpha_pi)|>attach()
#'epsilon_ba=.999*diag(4)+.001/3*(1-diag(4))
sampler_tau<-function(tau_vgb,pi_gs,epsilon_ba,n_vsa,
                      v=dim(tau_vgb)[1],
                      #s=dim(pi_gs)[2],
                      g=dim(pi_gs)[1],
                      g_neq_g=g_neq_g_f(g)){
  
  einsum::einsum(
    equation_string="vsgac,vsc->vga",
    log(
      einsum::einsum("v,gs,ac->vsgac",1:v,pi_gs,epsilon_ba)+
        einsum::einsum("a,gh,vhb,hs,bc->vsgac",rep(1,4),g_neq_g,tau_vgb,pi_gs,epsilon_ba)),
    n_vsa)|>
    plyr::aaply(1:2,function(x){exp(x-max(x))})|>
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
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha_pi=alpha_pi)|>attach()
#'epsilon_ba=.999*diag(4)+.001/3*(1-diag(4))
#'sampler_xi(tau_vgb,pi_gs,epsilon_ba)
sampler_xi<-function(tau_vgb,pi_gs,epsilon_ba,
                     v=dim(tau_vgb)[1],
                     s=dim(pi_gs)[2],
                     g=dim(pi_gs)[1]){
  ns<-  einsum::einsum("vsa,b,g->vsabg",n_vsa,rep(1,4),rep(1,g))
  ms<-einsum::einsum("vgb,gs,ba->vsabg",tau_vgb,pi_gs,epsilon_ba)
  nms=abind::abind(ns,ms,along=6)
  nm=nms[1,1,1,,,,drop=FALSE]|>abind::adrop(1:3)
  nms|>plyr::aaply(1:3,function(nm){
    rmultinom(1,size=nm[1,1,1],prob = nm[,,2])|>
      matrix(4)},.drop = FALSE)
  
}


sampler_xi_i<-function(i,tau_ivga,pi_igs,epsilon_iba,
                       v=dim(tau_vgb)[1],
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
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha_pi=alpha_pi)|>attach()
#'xi=sampler_xi(tau_vgb,pi_gs,epsilon_ba)
#'mu_vsag=mu_from_xi(xi)
#'rep_alpha_pi=rep(alpha_pi,g)
#'sampler_pi(mu_vsag,rep_alpha_pi)

sampler_pi<-function(mu_vsag,
                     rep_alpha_pi){
  
  einsum::einsum("vsag->sg",mu_vsag)|>
    plyr::aaply(1,function(x_g){
      x=dirmult::rdirichlet(1,alpha=rep_alpha_pi+x_g)
      x[x<.Machine$double.xmin]<-.Machine$double.xmin
      x/sum(x)})|>
    t()
}


sampler_pi_i<-function(mu_ivsag,
                       rep_alpha_pi){
  
  einsum::einsum("ivsag->isg",mu_ivsag)|>
    plyr::aaply(1:2,function(x_g){
      x=dirmult::rdirichlet(1,alpha=rep_alpha_pi+x_g)
      x[x<.Machine$double.xmin]<-.Machine$double.xmin
      x/sum(x)})|>
    aperm(c(1,3,2))
}

a_neq_b=g_neq_g_f(4)

#'@description sampler_epsilon_star
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha_pi=alpha_pi)|>attach()
#'
#'nu_vsab=nu_from_xi(xi)
#'alpha_epsilon=c(1,10)
#'
#'sampler_epsilon_star(nu_vsab,alpha_epsilon)

sampler_bar_epsilon_1_0<-function(alpha_epsilon){
  rbeta(1,shape1=alpha_epsilon[1],shape2=alpha_epsilon[2])}

sampler_bar_epsilon<-function(nu_vsab,alpha_epsilon){
  bar_epsilon1=rbeta(1,
                      shape1=alpha_epsilon[1]+einsum::einsum("vsab,ab->",nu_vsab,a_neq_b),
                      shape2=alpha_epsilon[2]+einsum::einsum("vsaa->",nu_vsab))}

sampler_epsilon_star<-function(nu_vsab,alpha_epsilon){
  bar_epsilon1=sampler_bar_epsilon(nu_vsab,alpha_epsilon)
  (1-4/3*bar_epsilon1)*diag(4)+(bar_epsilon1/3)}



#'@description SMC kernel
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha_pi=alpha_pi)|>attach()




desman_kernel<-function(n_vsa,tau_vgb,pi_gs,epsilon_ba,rep_alpha_pi,delta){
  xi_vsabg=sampler_xi(tau_vgb=tau_vgb,
                      pi_gs=pi_gs,
                      epsilon_ba = epsilon_ba)
  nu_vsab<-nu_from_xi(xi=xi_vsabg)
  mu_vsab<-mu_from_xi(xi=xi_vsabg)
  pi_gs<-sampler_pi(mu_vsag = mu_vsab,alpha_pi = rep_alpha_pi)
  w1<-smc_b_prime(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
  tau_vgb<-sampler_tau(n_vsa = n_vsa,tau_vgb = tau_vgb,pi_gs = pi_gs,epsilon_ba = epsilon_ba)
  w2<-smc_b_prime(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
  epsilon_ba<-sampler_epsilon_star(nu_vsab,delta)
  w3<-smc_b_prime(n_vsa_lambda,tau_ivga,pi_igs,epsilon_iba)
  list(tau_vgb=tau_vgb,pi_gs=pi_gs,epsilon_ba=epsilon_ba,w1=w1,w2=w2,w3=w3)
}

smc_kernel_i<-function(i,g,n_vsa,tau_ivga,pi_igs,epsilon_iba,rep_alpha_pi,delta,g_neq_g =g_neq_g_f(g)){
  xi_ivsabg=sampler_xi_i(i =i, tau_ivga=tau_ivga,
                         pi_igs=pi_igs,
                         epsilon_iba = epsilon_iba)
  nu_ivsab=nu_i_from_xi_i(xi=xi_ivsabg)
  mu_ivsab=mu_i_from_xi_i(xi=xi_ivsabg)
  pi_igs=sampler_pi_i(mu_ivsag = mu_ivsab,rep_alpha_pi = rep_alpha_pi)
  tau_ivga=sampler_tau_i(n_vsa = n_vsa,tau_ivga = tau_ivga,pi_igs = pi_igs,epsilon_iba = epsilon_iba,g_neq_g=g_neq_g )
  list(tau_ivga=tau_ivga,pi_igs=pi_igs,epsilon_iba=epsilon_iba)
}


