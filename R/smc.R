
g_neq_g_f=function(g){1-diag(g)}


#'@description sampler_tau
#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'sim_tau_pi_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
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
#'sim_tau_pi_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
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
#'sim_tau_pi_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
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
#'sim_tau_pi_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
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
#'sim_tau_pi_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()
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
#'@description SMC algorithm
#'@examples
#'n=1000;v=20;g=5;s=3;alpha0=.1
#'sim_tau_pi_n(v=v, g=g, s=s, n=n, error_rate = .001, alpha0=alpha0)|>attach()


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
