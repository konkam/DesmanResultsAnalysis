
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

#' @examples
#' sim_pi_gs(g = 3, s=5,alpha0=.1)
sim_pi_gs <- function(g, s, alpha0 = 1) {
  rdirichlet(alpha = rep(alpha0, g), n_samples = s) |> t()
}

epsilon_ba_f <- function(error_rate){diag(x = 1 - error_rate, nrow = 4) + 
  error_rate / 3 * (matrix(data = 1, nrow = 4, ncol = 4) - diag(x = 1, nrow = 4))}



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
  epsilon_ba <- epsilon_ba_f(error_rate)
  p_vsabg <- einsum::einsum(
    "vgb,gs,ba->vsabg",
    tau_vga,
    pi_gs,
    epsilon_ba
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
#' g = 5;v = 50; s = 3; n = 100; alpha0 = .1
#' sim_tau_pi_epsilon_n(g = 5, v = 50, s = 3, n = 100, alpha0 = 1)
sim_tau_pi_epsilon_n <- function(v, g, s, n, error_rate = .001, alpha0 = 1) {
  tau_vga <- sim_tau_vga(v = v, g = g)
  pi_gs <- sim_pi_gs(g = g, s = s, alpha0 = alpha0)
  epsilon_ba=epsilon_ba_f(error_rate)
  n_vsa <- sim_n_vsa(n = n, tau_vga = tau_vga, pi_gs = pi_gs, error_rate = error_rate)
  list(tau_vga = tau_vga, pi_gs = pi_gs, epsilon_ba=epsilon_ba,n_vsa = n_vsa)
}