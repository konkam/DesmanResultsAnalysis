nucleotides <- c("a", "c", "g", "t") |>
  (function(x) {
    x |> setNames(x)
  })()

# Helper functions
jags_sample_to_summary_tibble <- function(jags_sample) {
  jags_sample |>
    summary() |>
    (function(ss) {
      tibble(parname = rownames(ss)) |>
        bind_cols(as_tibble(ss))
    })()
}
## Data simulation

# Helper functions
#' @examples
# alpha_values <- c(2, 3, 4) # Parameters for the Dirichlet distribution
#' sample_dirichlet <- rdirichlet_onesmp(alpha_values)
#' print(sample_dirichlet)
rdirichlet_onesmp <- function(alpha) {
  # Generate K independent random samples from Gamma distributions
  y <- rgamma(length(alpha), shape = alpha, rate = 1)

  # Normalize the samples to obtain the Dirichlet variables
  x <- y / sum(y)

  return(x)
}
#' @examples
#' alpha_values <- c(2, 3, 4) # Parameters for the Dirichlet distribution
#' num_samples <- 5 # Number of samples
#' sample_dirichlet <- rdirichlet(alpha_values, n_samples = num_samples)
#' print(sample_dirichlet)
rdirichlet <- function(alpha, n_samples = 1) {
  # Generate n_samples independent random samples from Gamma distributions
  y <- matrix(rgamma(n_samples * length(alpha), shape = alpha, rate = 1), 
              ncol = length(alpha))

  # Normalize the samples to obtain the Dirichlet variables
  x <- t(apply(y, 1, function(row) row / sum(row)))

  return(x)
}


#' @examples
#' dna_string <- "AC"
#' result <- translate_dna_to_binary(dna_string)
#' print(result)
translate_dna_to_binary <- function(dna_string) {
  # Create a matrix to store the binary representation
  binary_matrix <- matrix(0,
    nrow = nchar(dna_string),
    ncol = 4,
    dimnames = list(NULL, nucleotides)
  )

  # Fill in the matrix based on the nucleotides in the DNA string
  for (i in 1:nchar(dna_string)) {
    nucleotide <- substr(dna_string, i, i)
    binary_matrix[i, nucleotide] <- 1
  }

  return(binary_matrix)
}


#' @examples
#' variants_string_matrix <- sim_variants_string_matrix(3, 12)
#' translate_dna_matrix_to_binary_array(variants)
translate_dna_matrix_to_binary_array <- function(variants_string_matrix) {
  variants_string_matrix |>
    tolower() |>
    plyr::aaply(1:2, `==`, nucleotides) |>
    (`*`)(1)
}

#' @examples
#' variants_string_vector <- sim_variants_string_matrix(g = 3, v = 12) |> plyr::aaply(2, paste0, collapse = "")
#' translate_dna_string_vector_to_string_matrix(variants_string_vector)
translate_dna_string_vector_to_string_matrix <- function(variants_string_vector) {
  variants_string_vector |>
    tolower() |>
    plyr::aaply(1, function(x) {
      strsplit(x, "") |>
        unname() |>
        unlist()
    }) |>
    t() |>
    (function(x) {
      x |> (`dimnames<-`)(list(v = 1:dim(x)[1], g = 1:dim(x)[2]))
    })()
}


#' @examples
#' sim_variants_string_matrix(v = 12, g = 3)
sim_variants_string_matrix <- function(v, g) {
  plyr::raply(g, sample(nucleotides, size = v, replace = TRUE)) |>
    t() |>
    (`dimnames<-`)(list(v = 1:v, g = 1:g))
}
#' @examples
#' sim_tau_vga(g = 3, v = 12)
sim_tau_vga <- function(v, g) {
  sim_variants_string_matrix(v, g) |>
    translate_dna_matrix_to_binary_array() |>
    (`*`)(1)
}

sim_pi_gs <- function(g, s, alpha0 = 1) {
  rdirichlet(alpha = rep(alpha0, g), n_samples = s) |> t()
}

#' @params
#' @examples
#' sim_n_vsa(s = 3, n = 1000, v = 50, g = 5)
sim_n_vsa <- function(n,
                      tau_vga = NULL,
                      pi_gs = NULL,
                      v = NULL,
                      g = NULL,
                      s = NULL,
                      error_rate = .001,
                      alpha0 = 1) {
  if (is.null(tau_vga)) {
    tau_vga <- sim_variants(v = v, g = g)
  } else {
    v <- dim(tau_vga)[1]
    g <- dim(tau_vga)[2]
  }
  if (is.null(pi_gs)) {
    pi_gs <- sim_pi_gs(g = g, s = s, alpha0 = alpha0)
  } else {
    s <- dim(pi_gs)[2]
  }
  # Mean coverage is 20, which is pretty favourable
  n_vs <- rpois(n = v * s, lambda = 20) |> matrix(c(v, s)) 
  epsilon <- diag(x = 1 - error_rate, nrow = 4) + 
    error_rate / 3 * (matrix(data = 1, nrow = 4, ncol = 4) - diag(x = 1, nrow = 4))
  
  p_vsabg <- einsum::einsum(
    "vgb,gs,ba->vsabg",
    tau_vga,
    pi_gs,
    epsilon
  )

  p_vsa <- plyr::aaply(p_vsabg, 1:3, sum)

  np_vsa <- abind::abind(n_vs, p_vsa, along = 3)

  n_vsa <- plyr::aaply(np_vsa, 1:2, function(x) {
    rmultinom(n = 1, size = x[1], prob = x[2:5]) |> setNames(nucleotides)
  })

  n_vsa |>
    (function(x) {
      x |> (`dimnames<-`)(dim(x) |> lapply(seq_len) |> setNames(c("v", "s", "a")))
    })()
}
#' @examples
#' sim_tau_pi_n(g = 5, v = 50, s = 3, n = n, alpha0 = 1)
sim_tau_pi_n <- function(v, g, s, n, error_rate = .001, alpha0 = 1) {
  tau_vga <- sim_tau_vga(v = v, g = g)
  pi_gs <- sim_pi_gs(g = g, s = s, alpha0 = alpha0)
  n_vsa <- sim_n_vsa(n = n, tau_vga = tau_vga, pi_gs = pi_gs, error_rate = error_rate)
  list(tau_vga = tau_vga, pi_gs = pi_gs, n_vsa = n_vsa)
}

## Inference
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
  return(c(a, b) |> setNames(c("a", "b")))
}



### Model


#' @description
#' Creates model string for jags
#' @param tildeepsilon a numerical value. if NA, then the model uses a dirichlet prior.
#' @examples
#' tau_pi_n <- sim_tau_pi_n(v = 50, g = 5, s = 3, n = 1000, alpha0 = 1)
#' cat(model_string_f(NA))
model_string_f <- function(tildeepsilon = NA,
                                gs="jags") {
  if(gs==jags){
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
      "tildeepsilon~ddirch(c(aa, bb))
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
  real aa;        // Dirichlet prior parameter for epsilon (match)
  real bb;        // Dirichlet prior parameter for epsilon (mismatch)"},
"
}

parameters {
  simplex[G] pi_gs[S];                     // Mixing proportions for each sample and group
  simplex[2] tildeepsilon;                 // Base probabilities for matching/mismatching
  real<lower=0, upper=1> tau_vga[V, G, 4]; // Proportion of each variant that matches group and allele
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
  "tildeepsilon ~ dirichlet([aa, bb]); // Dirichlet prior on base match/mismatch probabilities"
},"
  // Likelihood
  for (v in 1:V) {
    for (s in 1:S) {
      nvs[v, s] ~ multinomial(p_vsa[v, s]);
    }
  }
}")}}



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
#' desman_fixed_variants(
#'   n_vsa = tau_pi_n$n_vsa,
#'   tau_vga = tau_pi_n$tau_vga,
#'   G = 12
#' )
desman_fixed_variants <- function(n_vsa,
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


  model_string <- model_string_f(tildeepsilon = tildeepsilon,gs=gs)
  data_list <- list(
    V = V,G = G,S = S,      n_vsa = n_vsa,
    tau_vga = tau_vga,
    alpha = rep(alpha0, G),
    #epsilon = 2 * tildeepsilon * diag(4) + (1 - tildeepsilon),
    nvs = n_vsa |> apply(MARGIN = c(1, 2), FUN = sum),
  )
  
  if (!is.na(tildeepsilon)) {
    data_list=c(data_list,list(tildeepsilon=tildeepsilon))
  }else{
    error_matrix_prior <- 
      error_matrix_prior_specification(
        error_rate = error_rate, 
        prior_std = prior_std)
    data_list=c(data_list,
                list(
                  aa = error_matrix_prior["a"],
                  bb = error_matrix_prior["b"]))}
  
  monitor= ic("pi_gs", f (is.na(tildeepsilon)) {"tildeepsilon"})
  
  # Compiling and producing posterior samples from the model.
  gibbs_samples <- if(gs="jags"){
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
