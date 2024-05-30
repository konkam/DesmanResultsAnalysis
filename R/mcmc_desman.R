#' @description
#' Creates model string for jags or stan
#' @param tildeepsilon a numerical value. if NA, then the model uses a Dirichlet prior.
#' @param gs a character string. If "jags", the gibbs sampler used will be jags, stan if "stan".
#' @examples
#' cat(mcmc_model_string_desman_f(gs="jags",  tildeepsilon=NA))
#' cat(mcmc_model_string_desman_f(gs="jags",  tildeepsilon=.1))
#' cat(mcmc_model_string_desman_f(gs="nimble",tildeepsilon=NA))
#' cat(mcmc_model_string_desman_f(gs="nimble",tildeepsilon=.1))
#' cat(mcmc_model_string_desman_f(gs="jags",  tildeepsilon=NA,block_tau=FALSE))
#' cat(mcmc_model_string_desman_f(gs="jags",  tildeepsilon=.1,block_tau=FALSE))
#' cat(mcmc_model_string_desman_f(gs="nimble",tildeepsilon=NA,block_tau=FALSE))
#' cat(mcmc_model_string_desman_f(gs="nimble",tildeepsilon=.1,block_tau=FALSE))
mcmc_model_string_desman_f <- function(gs="jags",
                                       tildeepsilon=NA,
                                       block_tau=TRUE) {
  paste0(
    if(gs=="jags"){
"model "},
if(is.element(gs,c("jags","nimble"))){paste0(
  "{
  # Likelihood
  for (v in 1:V){
    for (s in 1:S){
      n_vsa[v,s,1:4] ~ dmulti(p_vsa[v,s,1:4],nvs[v,s])
    }
  }
  
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
            tau_vga[v,g,a] =ifelse(tau_vg[v,g]==a,1,0)
          }"},
        if(gs=="nimble"){"
          tau_vga[v,g,1:4]~dmulti(prob=p_tau[1:4],size=1)"},"}")},
      if(block_tau){
          "
          tau_v[v] ~ dcat(p_tau[1:G4])
          for (g in 1:G){
            tau_vg[v,g]=1+trunc((tau_v[v]-1)/(4^(g-1)))-4*trunc((tau_v[v]-1)/(4^(g)))
            for (a in 1:4){
              tau_vga[v,g,a] =ifelse(tau_vg[v,g]==a,1,0)
            }}"},
  "
          for (g in 1:G){
            for (a in 1:4){
              mixed_variants[v, g, a] = sum(tau_vga[v,g,1:4], epsilon[1:4,a])
            }
          }
      }

  # Latent multinomial observation probability
  for (v in 1:V){
    for (s in 1:S){
      for (g in 1:G){
        for (a in 1:4){
          p_g[v, s, g, a] = pi_gs[g, s] * mixed_variants[v, g, a]
        }
      }
      for (a in 1:4){
        p_vsa[v, s, a] = sum(p_g[v, s,1:G, a]) # Sum over variants
      }
    }
  }

  # Prior
  for(g in 1:G){
  repalpha_pi[g]<-alpha_pi
  }
  for (s in 1:S){
    pi_gs[1:G, s] ~ ddirch(repalpha_pi[1:G])
  }",
    if (is.na(tildeepsilon)) {
      "
      tildeepsilon~ddirch(shape_epsilon)
      "
    },
    "
     for (a in 1:4){
      for (b in 1:4){
        epsilon[a,b] = (a!=b)*tildeepsilon[2]/3 +(a==b)*tildeepsilon[1]
        }
    }
}
")})}


#' @description
#' Run jags.
#' @param n_vsa an array of counts.
#' @param tau_vga a collection of variants
#' @param G an integer. If G is smaller than the number of variants in the variant bin, then the algorithm is ran on the minimum of G and the number of variants in the bin.
#' @param tildeepsilon a numerical value. if NA, then the model uses a dirichlet prior.
#' @param error_rate = 0.001 controls the dirichlet prior on tildeepsilon
#' @param prior_std = 0.01 controls the dirichlet prior on tildeepsilon
#' @examples
#' tau_pi_n <- sim_tau_pi_epsilon_n(v = 50, g = 5, s = 3, n = 1000, alpha_pi = 1)
#' gs="jags"
#' block_tau=FALSE
#' GS=mcmc_unrelaxed_run(nburnin = 10,niter = 50,nchains = 1, gs="nimble", n_vsa = tau_pi_n$n_vsa,   G = 5,block_tau=FALSE)
#' GS=mcmc_unrelaxed_run(nburnin = 10,niter = 50,nchains = 1, gs="nimble", n_vsa = tau_pi_n$n_vsa,   G = 5,block_tau=FALSE,tildeepsilon=.1)
#' GS=mcmc_unrelaxed_run(nburnin = 10,niter = 50,nchains = 1, gs="nimble", n_vsa = tau_pi_n$n_vsa,   G = 5,block_tau=TRUE)
#' GS=mcmc_unrelaxed_run(nburnin = 10,niter = 50,nchains = 1, gs="nimble", n_vsa = tau_pi_n$n_vsa,   G = 5,block_tau=TRUE ,tildeepsilon=.1)
#' GS=mcmc_unrelaxed_run(burnin = 10,sample = 50,n.chains = 1,gs="jags",   n_vsa = tau_pi_n$n_vsa,   G = 5,block_tau=FALSE)
#' GS=mcmc_unrelaxed_run(burnin = 10,sample = 50,n.chains = 1,gs="jags",   n_vsa = tau_pi_n$n_vsa,   G = 5,block_tau=FALSE,tildeepsilon=.1)
#' GS=mcmc_unrelaxed_run(burnin = 10,sample = 50,n.chains = 1,gs="jags",   n_vsa = tau_pi_n$n_vsa,   G = 5,block_tau=TRUE)
#' GS=mcmc_unrelaxed_run(burnin = 10,sample = 50,n.chains = 1,gs="jags",   n_vsa = tau_pi_n$n_vsa,   G = 5,block_tau=TRUE ,tildeepsilon=.1)
#' GS$mcmc[[1]]->X;
#' X[10000,paste0("pi_gs[",c(outer(1:12,1:3,paste,sep=",")),"]")]
#' A=matrix(NA,12,3) ;for (i in 1:12){for(j in 1:3){A[i,j]=X[10000,paste0("pi_gs[",i,",",j,"]")]}}
mcmc_unrelaxed_run <- function(n_vsa,
                                G,
                                gs="jags",
                                block_tau=FALSE,
                                tildeepsilon=.01,
                                shape_epsilon=c(100,1),
                                alpha_pi=1,
                                n_chains =2,
                                ...) {
  V <- dim(n_vsa)[1]
  S <- dim(n_vsa)[2]
  
  dimnames(n_vsa) <- lapply(dim(n_vsa), seq_len)
  
  
  model_string <- mcmc_model_string_desman_f(gs=gs,tildeepsilon=tildeepsilon,block_tau=block_tau)
  
  data_list <- list(
    n_vsa = n_vsa,
    nvs = n_vsa |> apply(MARGIN = c(1, 2), FUN = sum))
  
  
  constants <- list(
    V = V,
    S = S,
    G = G,
    alpha_pi = alpha_pi
  )
  tau_vga=sim_tau_vga(g = G,v = V)
  tau_vg=plyr::aaply(tau_vga,1:2,function(x){sum(x*(0:3))+1})
  tau_v=plyr::aaply(tau_vg,1,function(x){sum((x-1)*(4^(0:(G-1))))+1})
  
  inits <- list(pi_gs = matrix(1/G, nrow = G, ncol = S),
                tau_vga=tau_vga,
                tau_vg=tau_vg,
                tau_v=tau_v,
                p_vsa=array(1/4,c(V,S,4)))
  monitor= c("pi_gs","tau_vga","tildeepsilon")
  
  if(block_tau){
    constants=c(constants,list(G4=4^G))
  }
  
  
  if(is.na(tildeepsilon)){
    inits     <- c(inits,list(tildeepsilon=c(.99,.01)))
    constants <- c(constants,list(shape_epsilon = shape_epsilon))
    monitor   <-c(monitor,"tildeepsilon")
  }else{
    constants <- c(constants,list(tildeepsilon=c(1-tildeepsilon,tildeepsilon)))
  }
  
  
    
  
  # Compiling and producing posterior samples from the model.
  if(gs=="jags"){
    gibbs_samples <- 
    runjags::run.jags(
      model = model_string,
      data = c(data_list,constants),
      monitor = monitor,
      ...)
    }
  if(gs=="nimble"){
    nimble_code=eval(parse(text=paste0("nimble::nimbleCode(",model_string,")")))
    model <- nimble::nimbleModel(
      code = nimble_code, 
      data = data_list, 
      constants = constants, 
      inits = inits)
    compiled_model <- compileNimble(model)   
    mcmcConf <- configureMCMC(model)
    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(model)
    Cmcmc <- compileNimble(Rmcmc, project = model)
    
    gibbs_samples <- runMCMC(Cmcmc, ...)
      
    }
  
  
  
  
    if(gs=="stan"){
      rstan::stan(model_code = model_string, data = c(data_list,constants), ...)}
  gibbs_samples
}


mcmc_desman_run <- function(n_vsa,
                            G,
                            gs="jags",
                            tildeepsilon=.01,
                            n_chains =2,
                            ...) {
  mcmc_unrelaxed_run(n_vsa,
                                                  G=G,
                                                  gs="jags",
                                                  block_tau=FALSE,
                                                  tildeepsilon=tildeepsilon,
                                                  shape_epsilon=NULL,
                                                  alpha_pi=1,
                                                  n_chains =n_chains,
                                                  ...)}


mcmc_desman_run<-function(n_vsa,
                      seed=1,
                      i,
                      g,
                      alpha_pi=.1,
                      tildeepsilon=.01,
                      shape_epsilon=c(1,100),
                      init=NULL){
  smc_sampler(n_vsa,
                        seed=1,
                        n_plus=sum(n_vsa),
                        n_vsa_df=reorder_counts(n_vsa = n_vsa,seed=seed),
                        g,
                        i,
                        t_min=1,
                        t_max,
                        ess_min,
                        smc_kernel=desman_kernel,
                        tildeepsilon=.1,
                        shape_epsilon=c(1,1000),
                        alpha_pi=.1,
                        trace_all=TRUE,
                        .update_lambda=update_lambda,
                        init=NULL,
                        mcmc=FALSE)}
