### Model
#' @description
#' Creates model string for jags
#' @param tildeepsilon a numerical value. if NA, then the model uses a Dirichlet prior.
#' @param gs a character string. If "jags", the gibbs sampler used will be jags, stan if "stan."
#' @examples
#' cat(model_string_fixed_variants_f(NA))
model_string_fixed_variants_f <- 
  function(tildeepsilon = NA,
          gs="jags") {
  if(gs=="jags"){
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
  for (v in 1:V){
    for (g in 1:G){
      for (a in 1:4){
        mixed_variants[v, g, a] = inprod(tau_vga[v,g,], epsilon[,a])
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
        p_vsa[v, s, a] = sum(p_g[v, s, , a]) # Sum over variants
      }
    }
  }

  # Prior
  for (s in 1:S){
    pi_gs[1:G, s] ~ ddirch(alpha[1:G])
    }",
    if (is.na(tildeepsilon)) {
      "tildeepsilon~ddirch(shape_epsilon)
"
    },
    "
     for (a in 1:4){
      for (b in 1:4){
        epsilon[a,b] = (a!=b)*tildeepsilon[2]/3 +(a==b)*tildeepsilon[1]
        }
    }
}
"
  )
}else{paste0(
    "
data {
  int<lower=1> V; // Number of positions
  int<lower=1> S; // Number of samples
  int<lower=1> G; // Number of variants
  int nvs[V, S];  // Total count for each position-sample combination
  real alpha[G];  // Dirichlet prior parameters for pi_gs",
    if (!is.na(tildeepsilon)) {
      "
  real shape_epsilon[2];        // Dirichlet prior parameter for epsilon (match)"},
"
}

parameters {
  simplex[G] pi_gs[S];                     // Mixing proportions for each sample and group
  simplex[2] tildeepsilon;                 // Base probabilities for matching/mismatching
  int<lower=0, upper=1> tau_vga[V, G, 4]; // Inidicator of each variant that matches position and variant
}

transformed parameters {
  real epsilon[4, 4];           // Epsilon matrix with different probabilities for match/mismatch
  real mixed_variants[V, G, 4]; // Computed variants after mixing
  real p_g[V, S, G, 4];         // Probability grid for variants, samples, groups, alleles
  real p_vsa[V, S, 4];          // Probability grid for observed data
  
  // Construct epsilon based on tildeepsilon
  for (a in 1:4) {
    for (b in 1:4) {
      epsilon[a, b] = (a != b) * tildeepsilon[2] / 3 + (a == b) * tildeepsilon[1];
    }
  }
  
  // Compute mixed variants based on tau_vga and epsilon
  for (v in 1:V) {
    for (g in 1:G) {
      for (a in 1:4) {
        mixed_variants[v, g, a] = dot_product(tau_vga[v, g], epsilon[, a]);
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
    pi_gs[s] ~ dirichlet(alpha); // Dirichlet prior on mixing proportions
  }
  
    }",
if (is.na(tildeepsilon)) {
  "tildeepsilon ~ dirichlet(shape_epsilon); // Dirichlet prior on base match/mismatch probabilities"
},"
  // Likelihood
  for (v in 1:V) {
    for (s in 1:S) {
      nvs[v, s] ~ multinomial(p_vsa[v, s]);
    }
  }
}")}}


#### Prior specification for the error matrix

#' @description
#' Let's call \eqn{\tilde\epsilon} the error rate, which is of order 10e-3.
#' \deqn{\epsilon=(1-\tilde\epsilon)\times I_4+ \frac{\tilde\epsilon}3\times(J_4-I_4),  } where \eqn{(\tilde\epsilon\sim\mathrm{Beta}(a,b)} means that the prior expected value for the diagonal element \eqn{\tilde\epsilon=\epsilon_{1,1}} is \eqn{\frac{a}{a+b}} and its variance is \eqn{\frac{ab}{(a+b)^2(a+b+1)}}
#' Moment-matching specification for \eqn{a,b}:
#' \deqn{
#' \begin{align}
#' &\frac{a}{a+b} = 1 - \eta \implies b = a\frac{\eta}{1-\eta} \\
#' &\text{Var}(\epsilon_{1,1}) = \frac{ab}{(a+b)^2(a+b+1)} = \frac{(1-\eta)\eta}{\frac{a}{1-\eta}+1} = \frac{(1-\eta)^2\eta}{a+1-\eta}\implies a = \eta - 1 + \frac{(1-\eta)^2\eta}{\text{Var}(\epsilon_{1,1})}
#' \end{align}
#' }
error_matrix_prior_specification <- function(error_rate, prior_std) {
  a <- error_rate - 1 + (1 - error_rate)^2 * error_rate / prior_std^2
  b <- a * error_rate / (1 - error_rate)
  return(c(a, b))
}


#' @description
#' Run jags.
#' @param n_vsa an array of counts.
#' @param tau_vga a collection of variants
#' @param G an integer. If G is smaller than the number of variants in the variant bin, then the algorithm is ran on the minimum of G and the number of variants in the bin.
#' @param tildeepsilon a numerical value. if NA, then the model uses a dirichlet prior.
#' @param error_rate = 0.001 controls the dirichlet prior on tildeepsilon
#' @param prior_std = 0.01 controls the dirichlet prior on tildeepsilon
#' @examples
#' tau_pi_n <- sim_tau_pi_n(v = 50, g = 5, s = 3, n = 1000, alpha0 = 1)
#' mcmc_fixed_variants(
#'   n_vsa = tau_pi_n$n_vsa,
#'   tau_vga = tau_pi_n$tau_vga,
#'   G = 12
#' )
mcmc_fixed_variants <- function(n_vsa,
                                  tau_vga,
                                  G,
                                  gs="jags",
                                  tildeepsilon = NA,
                                  error_rate = 0.001,
                                  prior_std = 0.01,
                                  n_chains = n_chains,
                                  alpha0 = 1,
                                  ...) {
  V <- dim(n_vsa)[1]
  S <- dim(n_vsa)[2]
  G <- dim(tau_vga)[2]
  dimnames(n_vsa) <- lapply(dim(n_vsa), seq_len)


  model_string <- model_string_fixed_variants_f(tildeepsilon = tildeepsilon,gs=gs)
  data_list <- list(
    V = V,G = G,S = S,      n_vsa = n_vsa,
    tau_vga = tau_vga,
    alpha = rep(alpha0, G),
    #epsilon = 2 * tildeepsilon * diag(4) + (1 - tildeepsilon),
    nvs = n_vsa |> apply(MARGIN = c(1, 2), FUN = sum))
  
  if (!is.na(tildeepsilon)) {
    data_list=c(data_list,list(tildeepsilon=c(1-tildeepsilon,tildeepsilon)))
  }else{
    shape_epsilon <- 
      error_matrix_prior_specification(
        error_rate = error_rate, 
        prior_std = prior_std)
    data_list=c(data_list,
                list(shape_epsilon=shape_epsilon))}
  
  monitor= c("pi_gs", if (is.na(tildeepsilon)) {"tildeepsilon"})
  
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
