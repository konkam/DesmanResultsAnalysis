#'@description sampler_tau
#'@examples
#'g_neq_g_f(5)
g_neq_g_f=function(g){1-diag(g)}


#'@description sampler_tau
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, epsilon_bar_1 = .001, alpha_pi=alpha_pi)|>attach()
#'epsilon_ba=.999*diag(4)+.001/3*(1-diag(4))
#'sampler_tau(tau_vgb=tau_vgb,pi_gs,epsilon_ba,n_vsa,v=dim(tau_vgb)[1],g=dim(pi_gs)[1],g_neq_g=g_neq_g_f(g),alpha_tau=0,m_vbg=NULL)|>
#'translate_dna_binary_array_to_string_vector()
sampler_tau<-function(tau_vgb,
                      pi_gs,
                      epsilon_ba,
                      n_vsa,
                      v=dim(tau_vgb)[1],
                      g=dim(pi_gs)[1],
                      g_neq_g=g_neq_g_f(g),
                      alpha_tau=0,
                      m_vbg=NULL){
  if(alpha_tau==0){
  einsum::einsum(
    equation_string="vsgab,vsa->vgb",
    log(
      einsum::einsum("v,gs,ba->vsgab",rep(1,v),pi_gs,epsilon_ba)+
        einsum::einsum("b,gh,vhc,hs,ca->vsgab",rep(1,4),g_neq_g,tau_vgb,pi_gs,epsilon_ba)),
    n_vsa)|>
    plyr::aaply(1:2,function(x){exp(x-max(x))})|>
    plyr::aaply(1:2,function(xi){sample(nucleotides,1,prob = xi)})|>
    translate_dna_matrix_to_binary_array()
    }else{
      plyr::aaply(m_vbg,c(1,3),
                  function(x_b){dirmult::rdirichlet(1,alpha=rep_alpha_tau+x_b)})
    }
  
}


#'@description sampler_tau
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n_i(i=30,v=v, g=g, s=s, n=n, epsilon_bar_1 = .001, alpha_pi=alpha_pi)|>attach()
#'epsilon_ba=.999*diag(4)+.001/3*(1-diag(4))
#'sampler_tau_ivgb(tau_ivgb=tau_ivgb,pi_igs=pi_igs,epsilon_iba=epsilon_iba,n_vsa,v=dim(tau_vgb)[1],g=dim(pi_gs)[1],g_neq_g=g_neq_g_f(g),alpha_tau=0,m_ivbg=NULL)

sampler_tau_ivgb<-function(tau_ivgb,
                        pi_igs,
                        epsilon_iba,
                        n_vsa,
                        v=dim(tau_ivgb)[2],
                        #s=dim(pi_gs)[2],
                        g=dim(pi_igs)[2],
                        g_neq_g=g_neq_g_f(g),
                        alpha_tau=0,
                        m_ivbg=NULL){
  if(alpha_tau==0){
  exp(
    einsum::einsum(
      equation_string="ivsgab,vsa->ivgb",
      log(
        einsum::einsum("v,igs,iba->ivsgab",rep(1,v),pi_igs,epsilon_iba)+
          einsum::einsum("b,gh,ivhc,ihs,ica->ivsgab",rep(1,4),g_neq_g,tau_ivgb,pi_igs,epsilon_iba)),
      n_vsa))|>
    plyr::aaply(1:3,function(xi){sample(nucleotides,1,prob = xi)})|>
    plyr::aaply(1,translate_dna_matrix_to_binary_array)
  }else{
    plyr::aaply(m_ivbg,c(1:2,4),
                function(x_b){dirmult::rdirichlet(1,alpha=rep_alpha_tau+x_b)})
  }
}




#'@description sampler_tau
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n_i(i=30,v=v, g=g, s=s, n=n, epsilon_bar_1 = .001, alpha_pi=alpha_pi)|>attach()
#'epsilon_ba=.999*diag(4)+.001/3*(1-diag(4))
#'sampler_tilde_tau_ivgb(rep_alpha_tau=rep(.1,4),m_ivba=)

sampler_tilde_tau_ivgb<-function(rep_alpha_tau,
                           m_ivgb){
    plyr::aaply(m_ivgb,c(1:3),
                function(x_b){dirmult::rdirichlet(1,alpha=rep_alpha_tau+x_b)})|>
    namedims("ivgb")}


sampler_tilde_rho_ivga<-function(rep_alpha_rho,
                                m_ivga){
  plyr::aaply(m_ivga,c(1:3),
              function(x_a){dirmult::rdirichlet(1,alpha=rep_alpha_rho+x_a)})|>
    namedims("ivga")}




#'@description sampler_mu_nu
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, epsilon_bar_1 = .001, alpha_pi=alpha_pi)|>attach()
#'epsilon_ba=.999*diag(4)+.001/3*(1-diag(4))
#'sampler_m(tau_vgb,pi_gs,epsilon_ba)
sampler_m<-function(tau_vgb,pi_gs,epsilon_ba,
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

#'@examples
#'index="ivsabg"
#'x=array(1:4^6,dim=rep(4,6))
namedims<-function(x,index){
  dimnames(x)<-plyr::llply(dim(x),seq_len)
  names(dimnames(x))<-strsplit(index, split = "")[[1]]
  which=is.element(names(dimnames(x)),c("a","b"))
  dimnames(x)[which]<-rep(list(nucleotides),sum(which));
  x}
#'@description sampler_mu_nu
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'theta=sim_tau_pi_epsilon_n_i(i=4,v=v, g=g, s=s, n=n, epsilon_bar_1 = .001, alpha_pi=alpha_pi)
#'obs=sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, epsilon_bar_1 = .001, alpha_pi=alpha_pi)
#'sampler_m_ivsabg(theta$chi_ivsabg,obs$n_vsa)

sampler_m_ivsabg<-function(chi_ivsabg,
                           n_vsa,
                           i=dim(chi_ivsabg)[1],
                           v=dim(chi_ivsabg)[2],
                           s=dim(chi_ivsabg)[3],
                           g=dim(chi_ivsabg)[6],
                           grid_ivsa=expand.grid(i=1:i,v=1:v,s=1:s,a=1:4)){
  
  
  
  plyr::maply(grid_ivsa,function(i,v,s,a){
    rmultinom(1,size=n_vsa[v,s,a],
              prob = chi_ivsabg[i,v,s,a,,])})|>
    array(dim=dim(chi_ivsabg))|>
    namedims("ivsabg")
  }


#'@examples
#'i=4;v=6; g=5; s=3; n=50;
#'sim_i=sim_tau_pi_epsilon_n_i(i=i,v=v, g=g, s=s, n=50, epsilon_bar_1 = .001, alpha_pi = 1)
#'chi_ivsag=sim_i$chi_ivsag
#'sim=sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=50)
#'n_vsa=sim$n_vsa
#'sampler_m_ivsag(chi_ivsag=chi_ivsag,n_vsa=n_vsa)
sampler_m_ivsag<-function(chi_ivsag,
                      n_vsa,
                      i=dim(chi_ivsag)[1],
                      v=dim(chi_ivsag)[2],
                      s=dim(chi_ivsag)[3],
                      g=dim(chi_ivsag)[5],
                      grid_ivsa=expand.grid(i=1:i,v=1:v,s=1:s,a=1:4)){
  plyr::maply(grid_ivsa,function(i,v,s,a){
    rmultinom(1,size=n_vsa[v,s,a],
              prob = chi_ivsag[i,v,s,a,])})|>
    namedims("ivsag")
  
}
m_vsab_from_m<-function(xi){einsum::einsum("vsabg->vsab",xi)}
m_vsag_from_m<-function(xi){einsum::einsum("vsabg->vsag",xi)}

m_ivsab_from_m_i<-function(xi){einsum::einsum("ivsabg->ivsab",xi)}
m_ivsag_from_m<-function(xi){einsum::einsum("ivsabg->ivsag",xi)}

#'@description sampler_pi
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, epsilon_bar_1 = .001, alpha_pi=alpha_pi)|>attach()
#'xi=sampler_m(tau_vgb,pi_gs,epsilon_ba)
#'mu_vsag=m_vsag_from_m(xi)
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


sampler_pi_igs<-function(m_igs,
                       rep_alpha_pi){
  plyr::aaply(m_igs,c(1,3),
              function(x_g){dirmult::rdirichlet(1,alpha=rep_alpha_pi+x_g)})|>
    aperm(c(1,3,2))|>
    namedims("igs")
}

a_neq_b=g_neq_g_f(4)

#'@description sampler_epsilon_star
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, epsilon_bar_1 = .001, alpha_pi=alpha_pi)|>attach()
#'
#'nu_vsab=m_vsab_from_m(xi)
#'alpha_epsilon=c(1,10)
#'
#'sampler_bar_epsilon_1_0(nu_vsab,alpha_epsilon)
#'sampler_bar_epsilon(nu_vsab,alpha_epsilon)

sampler_epsilon_iba<-function(m_iba,rep_alpha_epsilon){
  plyr::aaply(m_iba,1:2,
              function(x_a){dirmult::rdirichlet(1,alpha=rep_alpha_epsilon+x_a)})|>
    namedims("iba")
}


sampler_bar_epsilon_i1<-function(m_ianeqb,m_iaa,alpha_bar_epsilon_1){
  plyr::aaply(1:nrow(m_ianeqb),1,
              function(i){
  rbeta(1,
                     shape1=alpha_bar_epsilon_1[1]+m_ianeqb,
                     shape2=alpha_bar_epsilon_1[2]+m_iaa)})}


sampler_bar_epsilon_1_0<-function(alpha_epsilon){
  rbeta(1,shape1=alpha_epsilon[1],shape2=alpha_epsilon[2])}

sampler_bar_epsilon<-function(nu_vsab,alpha_epsilon){
  bar_epsilon1=rbeta(1,
                      shape1=alpha_epsilon[1]+einsum::einsum("vsab,ab->",nu_vsab,a_neq_b),
                      shape2=alpha_epsilon[2]+einsum::einsum("vsaa->",nu_vsab))}





sampler_epsilon_star<-function(nu_vsab,alpha_epsilon){
  bar_epsilon1=sampler_bar_epsilon(nu_vsab,alpha_epsilon)
  (1-4/3*bar_epsilon1)*diag(4)+(bar_epsilon1/3)}

