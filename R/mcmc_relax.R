
#' @description
#' Creates model string for jags or stan
#' @param gs a character string. If "jags", the gibbs sampler used will be jags, stan if "stan."
#' @param block_rho a boolean. If TRUE, dictionnary is sampled simultaneously for all variants at each position.
#' @examples
#' cat(mcmc_relax_model_string_fixed_variants_f(block_rho=TRUE))
#' 
mcmc_relax_model_string_fixed_variants_f <- function(block_rho=TRUE,gs="jags") {
  paste0(
    "
model {
  # Likelihood
  for (v in 1:V){
    for (s in 1:S){
      n_vsa[v,s,] ~ dmulti(p_vsa[v,s,], nvs[v,s])
    }
  }
  
  # Mangled variants
  for (v in 1:V){",
    if(block_rho){"block {
        "},
    "for (g in 1:G){
        tau[v, g,1:4 ] ~ ddirich(rep(alpha_rho,4))
        rho[v,g,1:4]=(1-4*tildeepsilon/3)*tau[v, g, 1:4]+(tildeepsilon/3)*rep(1,4)
    }",
    if(block_rho){"
        }"},
    "}
  }
  # Latent multinomial observation probability
  for (v in 1:V){
    for (s in 1:S){
      for (g in 1:G){
        for (a in 1:4){
          p_g[v, s, g, a] = pi_gs[g, s] * rho[v, g, a]
        }
      }
      for (a in 1:4){
        p_vsa[v, s, a] = sum(p_g[v, s, , a]) # Sum over variants
      }
    }
  }

  # Prior
  for (s in 1:S){
    pi_gs[1:G, s] ~ ddirch(rep(alpha_pi,G))
  }
  
  tildeepsilon~ddirch(shape_epsilon)")}


#' @description
#' Run jags.
#' @param n_vsa an array of counts.
#' @param tau_vga a collection of variants
#' @param G an integer, maximum number of variants.
#' @param shape_rho a length two vector.
#' @param alpha_pi
#' @examples
#' tau_pi_n <- sim_tau_pi_n(v = 50, g = 5, s = 3, n = 1000, alpha0 = 1)
#' gs="jags"
#' block_tau=FALSE
#' shape_rho=c(1,100)
#' GS=mcmc_very_relax_run(
#'   n_vsa = tau_pi_n$n_vsa,
#'   G = 12
#' )
#' GS$mcmc[[1]]->X;
#' X[10000,paste0("pi_gs[",c(outer(1:12,1:3,paste,sep=",")),"]")]
#' A=matrix(NA,12,3) ;for (i in 1:12){for(j in 1:3){A[i,j]=X[10000,paste0("pi_gs[",i,",",j,"]")]}}
mcmc_relax_run <- function(n_vsa,
                                G,
                                gs="jags",
                                shape_rho=c(1,100),
                                block_rho=FALSE,
                                alpha_pi=1,
                                n_chains =2,
                                ...) {
  V <- dim(n_vsa)[1]
  S <- dim(n_vsa)[2]
  
  dimnames(n_vsa) <- lapply(dim(n_vsa), seq_len)
  
  
  model_string <- mcmc_very_relax_model_string_fixed_variants_f(gs=gs,block_rho = block_rho)
  data_list <- list(
    V = V,G = G,S = S,      
    n_vsa = n_vsa,
    shape_rho=shape_rho,alpha_pi=alpha_pi,
    #epsilon = 2 * tildeepsilon * diag(4) + (1 - tildeepsilon),
    nvs = n_vsa |> apply(MARGIN = c(1, 2), FUN = sum))
  
  
  monitor= c("pi_gs","tau","rho")
  
  # Compiling and producing posterior samples from the model.
  gibbs_samples <- if(gs=="jags"){
    runjags::run.jags(
      model = model_string,
      data = data_list,
      monitor = monitor,
      ...
    )}else{
      rstan::stan(model_code = model_string, data = data_list, ...)
    }
  gibbs_samples
}