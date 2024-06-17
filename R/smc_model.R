#' @description
#' Creates model string for jags or stan
#' @param bar_epsilon a numerical value. if NA, then the model uses a Dirichlet prior.
#' @param gs a character string. If "jags", the gibbs sampler used will be jags, stan if "stan".
#' @examples
#' cat(model_string_f(gs="jags",fixed_bar_epsilon=FALSE))
model_string_f <- 
  function(gs="jags",
           fixed_bar_epsilon=FALSE,
           constrained_epsilon_matrix=TRUE,
           block_tau=TRUE,
           fixed_tau=FALSE,
           relax_tau=FALSE,
           relax_rho=FALSE,
           fixed_alpha_rho=FALSE) {
  paste0(
    if(gs=="jags"){
"model "},
  "{
  # Likelihood
  for (v in 1:V){
    for (s in 1:S){
      n_vsa[v,s,1:4] ~ dmulti(p_vsa[v,s,1:4],n_vs[v,s])
    }
  }",
  #---------------------------------------tau
  if((!relax_rho)){paste0(
  "
  # Mangled variants",
  if(!fixed_tau){
    if(relax_tau){
    "
    for(a in 1:4){
    rep_alpha_tau[a]=alpha_tau}
    for(v in 1:V){
      for(g in 1:G){
        tau_vgb[v,g,1:4]~ddirch(rep_alpha_tau[1:4])  
      }
    }  
    "
      }else{
    paste0(
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
            }}"},  "}")}},
  #---------------------------------------epsilon
  if (!fixed_bar_epsilon&constrained_epsilon_matrix){
    "
      bar_epsilon~ddirch(alpha_bar_epsilon)
      bar_epsilon_1=bar_epsilon[1]
      "
  },
  if (fixed_bar_epsilon|constrained_epsilon_matrix){
    "
  for (a in 1:4){
      for (b in 1:4){
        epsilon[a,b] = (a!=b)*bar_epsilon[2]/3 +(a==b)*bar_epsilon[1]
        }
    }
      "
  },
  if (!fixed_bar_epsilon&!constrained_epsilon_matrix){
    "
  for (a in 1:4){
        epsilon[a,1:4] = ddirch(alpha_epsilon)
    }
    "})},
#---------------------------------------rho
if(!relax_rho){
"for(v in 1:V){
      for (g in 1:G){
            for (a in 1:4){
              rho_vga[v, g, a] = inprod(tau_vgb[v,g,1:4], epsilon[1:4,a])
            }
          }
  }"},
if(relax_rho&!fixed_alpha_rho){
  "
    alpha_rho~dbeta(kappa_rho[1],kappa_rho[2])
    for(a in 1:4){rep_alpha_rho[a]=alpha_rho}"},
if(relax_rho){"
    for(v in 1:V){
      for (g in 1:G){
              rho_vga[v, g, 1:4]~ddirch(rep_alpha_rho) 
            
          }
  }
    "
},
#---------------------------------------n 
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
  }
  
  
  for (s in 1:S){
    pi_gs[1:G, s] ~ ddirch(rep_alpha_pi[1:G])
  }
}
")}

stan_model_string_generic_f <- function(gs="jags",
                                        fixed_bar_epsilon=FALSE,
                                        constrained_epsilon_matrix=TRUE,
                                        block_tau=TRUE,
                                        fixed_tau=FALSE,
                                        relax_tau=FALSE,
                                        relax_rho=FALSE,
                                        fixed_alpha_pi=FALSE){paste0(
                                          "
data {
  int<lower=1> V; // Number of positions
  int<lower=1> S; // Number of samples
  int<lower=1> G; // Number of variants
  int nvs[V, S];  // Total count for each position-sample combination
  real alpha[G];  // Dirichlet prior parameters for pi_gs",
                                          if (!is.null(bar_epsilon_1)) {
                                            "
  real shape_epsilon[2];        // Dirichlet prior parameter for epsilon (match)"},                                         "
}

parameters {
  simplex[G] pi_gs[S];                     // Mixing proportions for each sample and group
  simplex[2] bar_epsilon;                 // Base probabilities for matching/mismatching
  int<lower=0, upper=1> tau_vgb[V, G, 4]; // Inidicator of each variant that matches position and variant
}

transformed parameters {
  real epsilon[4, 4];           // Epsilon matrix with different probabilities for match/mismatch
  real mixed_variants[V, G, 4]; // Computed variants after mixing
  real p_g[V, S, G, 4];         // Probability grid for variants, samples, groups, alleles
  real p_vsa[V, S, 4];          // Probability grid for observed data
  
  // Construct epsilon based on bar_epsilon
  for (a in 1:4) {
    for (b in 1:4) {
      epsilon[a, b] = (a != b) * bar_epsilon[2] / 3 + (a == b) * bar_epsilon[1];
    }
  }
  
  // Compute mixed variants based on tau_vgb and epsilon
  for (v in 1:V) {
    for (g in 1:G) {
      for (a in 1:4) {
        mixed_variants[v, g, a] = dot_product(tau_vgb[v, g], epsilon[, a]);
      }
    }
  }
  
  // Calculate p_g and p_vsa using mixed_variants and pi_gs
  for (v in 1:V) {
    for (s in 1:S) {
      for (g in 1:G) {
        for (a in 1:4) {
          p_g[v, s, g, a] = pi_gs[s, g] * mixed_variants[v, g, a];
        }
      }
      for (a in 1:4) {
        p_vsa[v, s, a] = sum(p_g[v, s, , a]); // Sum over G
      }
    }
  }
}

model {
  // Prior distributions
  for (s in 1:S) {
    pi_gs[s] ~ ddirch(rep_alpha_pi); // Dirichlet prior on mixing proportions
  }
  
    }",
if (is.null(bar_epsilon_1)) {
  "bar_epsilon ~ ddirch(alpha_bar_epsilon); // Dirichlet prior on base match/mismatch probabilities
  bar_epsilon_1=bar_epsilon[1]"
},"
  // Likelihood
  for (v in 1:V) {
    for (s in 1:S) {
      nvs[v, s] ~ multinomial(p_vsa[v, s]);
    }
  }
}")} 

model_string_generic_f <-function(gs="jags",
                                  fixed_bar_epsilon=FALSE,
                                  constrained_epsilon_matrix=TRUE,
                                  block_tau=TRUE,
                                  fixed_tau=FALSE,
                                  relax_tau=FALSE,
                                  relax_rho=FALSE,
                                  fixed_alpha_pi=FALSE){
  
  if(gs=="stan"){
    stan_model_string_generic_f(
      fixed_bar_epsilon=fixed_bar_epsilon,
      constrained_epsilon_matrix=constrained_epsilon_matrix,
      block_tau=block_tau,
      fixed_tau=fixed_tau,
      relax_tau=relax_tau,
      relax_rho=relax_rho,
      fixed_alpha_pi=fixed_alpha_pi) 
  }else{
    model_string_f(
      fixed_bar_epsilon=fixed_bar_epsilon,
      constrained_epsilon_matrix=constrained_epsilon_matrix,
      block_tau=block_tau,
      fixed_tau=fixed_tau,
      relax_tau=relax_tau,
      relax_rho=relax_rho,
      fixed_alpha_pi=fixed_alpha_pi) 
  }
}