
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



## Inference
