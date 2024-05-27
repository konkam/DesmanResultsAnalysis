
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
#' dna_string <- "ActGG"
#' result <- translate_dna_to_binary(dna_string)
#' print(result)
translate_dna_to_binary <- function(dna_string) {
  # Create a matrix to store the binary representation
  binary_matrix <- matrix(0,
                          nrow = dna_string|>nchar(),
                          ncol = 4,
                          dimnames = list(NULL, nucleotides)
  )
  
  # Fill in the matrix based on the nucleotides in the DNA string
  for (i in 1:nchar(dna_string)) {
    nucleotide <- dna_string|>tolower()|>substr(i, i)
    binary_matrix[i, nucleotide] <- 1
  }
  
  return(binary_matrix)
}


#' @examples
#' variants_string_matrix <- sim_variants_string_matrix(3, 12)
#' translate_dna_matrix_to_binary_array(variants_string_matrix)
translate_dna_matrix_to_binary_array <- function(variants_string_matrix) {
  variants_string_matrix |>
    tolower() |>
    plyr::aaply(1:2, `==`, nucleotides) |>
    (`*`)(1)
}

#' @examples
#' variants_string_vector <- 
#' sim_variants_string_matrix(g = 3, v = 12) |> 
#' plyr::aaply(2, paste0, collapse = "")
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
#' variants_string_vector <- 
#' sim_variants_string_matrix(g = 3, v = 12) |> 
#' plyr::aaply(2, paste0, collapse = "")
#' translate_dna_string_vector_to_string_matrix(variants_string_vector)|>
#' translate_dna_matrix_to_binary_array()|>
#' translate_dna_binary_array_to_string_vector()
translate_dna_binary_array_to_string_vector <- 
  function(variants_binary_array) {
    variants_binary_array |>
      plyr::aaply(1:2,function(x){nucleotides[x==1]},.drop=FALSE)|>
      abind::adrop(3)|>
      plyr::aaply(2,paste,collapse="")
  }

## Inference
