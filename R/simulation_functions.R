
#' @examples
#' sim_variants_string_matrix(v = 12, g = 3)
nucleotides <- c("a", "c", "g", "t") |>
  (function(x) {
    x |> setNames(x)
  })()

#' @examples
#' sim_variants_string_matrix(v = 12, g = 3)
sim_variants_string_matrix <- function(v, g) {
  plyr::raply(g, sample(nucleotides|>unname(), size = v, replace = TRUE)) |>
    t() |>
    (`dimnames<-`)(list(v = 1:v, g = 1:g))
}
#' @examples
#' sim_tau_vgb(g = 3, v = 12)
sim_tau_vgb <- function(v, g) {
  sim_variants_string_matrix(v, g) |>
    translate_dna_matrix_to_binary_array() |>
    (`*`)(1)
}

#' @examples
#' sim_pi_gs(g = 3, s=5,alpha_pi=.1)
sim_pi_gs <- function(g, s, alpha_pi = 1) {
  rdirichlet(alpha = rep(alpha_pi, g), n_samples = s) |> 
    t()|>(`dimnames<-`)(list(g=1:g,s=1:s))
}

epsilon_ba_f <- function(error_rate){diag(x = 1 - error_rate, nrow = 4) + 
  error_rate / 3 * (matrix(data = 1, nrow = 4, ncol = 4) - diag(x = 1, nrow = 4))}



#' @params
#' @examples
#' sim_n_vsa(s = 3, n = 1000, v = 50, g = 5)
sim_n_vsa <- function(n=1000,# expeted sample size
                      tau_vgb = NULL,
                      pi_gs = NULL,
                      v = NULL,
                      g = NULL,
                      s = NULL,
                      error_rate = .001,
                      alpha_pi = 1) {
  if (is.null(tau_vgb)) {
    tau_vgb <- sim_variants(v = v, g = g)
  } else {
    v <- dim(tau_vgb)[1]
    g <- dim(tau_vgb)[2]
  }
  if (is.null(pi_gs)) {
    pi_gs <- sim_pi_gs(g = g, s = s, alpha_pi = alpha_pi)

  } else {
    s <- dim(pi_gs)[2]
  }
  
  # Mean coverage is 20, which is pretty favourable
  n_vs <- if(is.null(n)){rpois(n = v * s, lambda = n) |> array(c(v, s))}else{array(n,c(v,s))}
  epsilon_ba <- epsilon_ba_f(error_rate)
  p_vsabg <- einsum::einsum(equation_string = "vgb,gs,ba->vsabg",
    tau_vgb,
    pi_gs,
    epsilon_ba
  )
  
  p_vsa <- plyr::aaply(p_vsabg, 1:3, sum,.drop = FALSE)|>abind::adrop(drop=4)
  
  
  np_vsa <- abind::abind(n_vs, p_vsa, along = 3)
  
  n_vsa <- plyr::aaply(np_vsa, 1:2, function(x) {
    rmultinom(n = 1, size = x[1], prob = x[2:5]) |> setNames(nucleotides)
  },.drop=FALSE)|>abind::adrop(drop=4)
  
  n_vsa |>
    (function(x) {
      x |> (`dimnames<-`)(dim(x) |> lapply(seq_len) |> setNames(c("v", "s", "a")))
    })()
}

#' @examples
#' g = 5;v = 50; s = 3; n = 100; alpha_pi = .1
#' sim_tau_pi_epsilon_n(g = 5, v = 50, s = 3, n = 100, alpha_pi = 1)
sim_tau_pi_epsilon_n <- function(v, g, s, n, error_rate = .001, alpha_pi = 1) {
  tau_vgb <- sim_tau_vgb(v = v, g = g)
  pi_gs <- sim_pi_gs(g = g, s = s, alpha_pi = alpha_pi)
  reorder_g=order(pi_gs|>plyr::aaply(1,sum),decreasing = TRUE)|>
    (function(x){x|>setNames(x)})()
  pi_gs=pi_gs[reorder_g,]|>(`dimnames<-`)(list(g=1:g,s=1:s))
  tau_vgb=tau_vgb[,reorder_g,]|>(`dimnames<-`)(list(v=1:v,g=1:g,b=nucleotides))
  epsilon_ba=epsilon_ba_f(error_rate)
  n_vsa <- sim_n_vsa(n = n, tau_vgb = tau_vgb, pi_gs = pi_gs, error_rate = error_rate)
  list(tau_vgb = tau_vgb, pi_gs = pi_gs, epsilon_ba=epsilon_ba,n_vsa = n_vsa)
}
