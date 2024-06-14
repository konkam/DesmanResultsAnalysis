#' @description
#' Creates model string for jags or stan
#' @param bar_epsilon a numerical value. if NA, then the model uses a Dirichlet prior.
#' @param gs a character string. If "jags", the gibbs sampler used will be jags, stan if "stan".
#'model_string_failed <- 
  function(gs="jags",
           bar_epsilon=NA,
           block_tau=TRUE,
           latent_n=FALSE,
           fixed_tau=FALSE,
           alpha_tau=0.01,
           alpha_epsilon=NA,
           relax=FALSE,
           kappa_rho=NULL) {
    paste0(
      if(gs=="jags"){
        "model "},
      paste0(
        if(!latent_n){
          "{
  # Likelihood
  for (v in 1:V){
    for (s in 1:S){
      n_vsa[v,s,1:4] ~ dmulti(p_vsa[v,s,1:4],n_vs[v,s])
    }
  }"},
        if(!latent_n&is.null(kappa_rho)){
          "
  for (v in 1:V){
    for (s in 1:S){
      for(a in 1:4){
        m_vsa[v,s,a,1:G4] ~ dmulti(chi_vs[v,s,a,1:G4],n_vsa[v,s,a])
        n_vsa[v,s,a] ~ dmulti([v,s,1:4],n_vs[v,s])
      }
    }
  }
    "}
        ,
        if(is.null(kappa_rho)){paste0(
          "
  # Mangled variants",
          if(!block_tau){"
  for (a in 1:4){
    p_tau[a]<- 1/4}"},
          if(block_tau){"
  for (a in 1:G4){
    p_tau[a]<- 1/(G4)}"},
          "
  for (v in 1:V){",
          if(!block_tau){
            paste0("
    for (g in 1:G){",
                   if(gs=="jags"){"
          tau_vg[v,g] ~ dcat(p_tau)
          for (a in 1:4){
            tau_vgb[v,g,a] =ifelse(tau_vg[v,g]==a,1,0)
          }"},
                   if(gs=="nimble"){"
          tau_vgb[v,g,1:4]~dmulti(prob=p_tau[1:4],size=1)"},"}")},
          if(block_tau){
            "
          tau_v[v] ~ dcat(p_tau[1:G4])
          for (g in 1:G){
            tau_vg[v,g]=1+trunc((tau_v[v]-1)/(4^(g-1)))-4*trunc((tau_v[v]-1)/(4^(g)))
            for (a in 1:4){
              tau_vgb[v,g,a] =ifelse(tau_vg[v,g]==a,1,0)
            }}"},
          if(!latent_n){
            "
          for (g in 1:G){
            for (a in 1:4){
              rho_vga[v, g, a] = inprod(tau_vgb[v,g,1:4], epsilon[1:4,a])
            }
          }
  }"},
          if(latent_n){
            "
          for (v in 1:V){
            for (s in 1:S){
              for (a in 1:4){
                for (b in 1:4){
                  for (g in 1:G){
                chi_vs[v, s, ((a-1)*4*G)+(G*(b-1))+g] = tau_vgb[v,g,b]* epsilon[b,a]*pi_gs[g,s])
            }
          }}}
  }"},
          
          if (is.na(bar_epsilon_1)&(length(alpha_epsilon)==2)){
            "
      bar_epsilon~ddirch(alpha_epsilon)
      "
          },
          if (!is.na(bar_epsilon_1)|(length(alpha_epsilon)==2)){
            "
  for (a in 1:4){
      for (b in 1:4){
        epsilon[a,b] = (a!=b)*bar_epsilon[2]/3 +(a==b)*bar_epsilon[1]
        }
    }
      "
          },
          if (!is.na(bar_epsilon_1)|(length(alpha_epsilon)==2)){
            "
  for (a in 1:4){
        epsilon[a,1:4] = ddirch(alpha_epsilon)
    }
    "})},
        if(!is.null(kappa_rho)){
          "
    alpha_rho~dbeta(kappa_rho[1],kappa_rho[2])
    for(a in 1:4){rep_alpha_rho[a]=alpha_rho}
    alpha_rho~ddirch(rep_alpha_rho)
    "
        },
        if(!latent_n){
          "
  # Latent multinomial observation probability
  for (v in 1:V){
    for (s in 1:S){
      for (g in 1:G){
        for (a in 1:4){
          p_g[v, s, g, a] = pi_gs[g, s] * rho_vga[v, g, a]
        }
      }
      for (a in 1:4){
        p_vsa[v, s, a] = sum(p_g[v, s,1:G, a]) # Sum over variants
      }
    }
  }"},
        "
  # Prior
  for(g in 1:G){
  rep_alpha_pi[g]<-alpha_pi
  }
  for (s in 1:S){
    pi_gs[1:G, s] ~ ddirch(rep_alpha_pi[1:G])
  }
}
"))}