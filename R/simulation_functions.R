
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
#' sim_tau_ivgb(g = 3, v = 12,i=2)
sim_tau_ivgb <- function(i,v, g) {
  plyr::raply(i,sim_tau_vgb(v,g),.drop = FALSE)|>
    (function(x){dimnames(x)[1]<-list(i=1:i);names(dimnames(x))[1]<-"i";x})()
}

#' @examples
#' sim_pi_gs(g = 3, s=5,alpha_pi=.1)
sim_pi_gs <- function(g, s, alpha_pi = 1) {
  rdirichlet(alpha = rep(alpha_pi, g), n_samples = s) |> 
    t()|>(`dimnames<-`)(list(g=1:g,s=1:s))
}
#' @examples
#' sim_pi_igs(i=3,g = 3, s=5,alpha_pi=.1)
sim_pi_igs <- function(i,g, s, alpha_pi = 1) {
  plyr::raply(i,sim_pi_gs(g=g,s=s,alpha_pi = alpha_pi),.drop = FALSE)|>
    (function(x){dimnames(x)[1]<-list(i=1:i);names(dimnames(x))[1]<-"i";x})()
}
#' @examples
#' epsilon_ba_f(.1)
epsilon_ba_f <- function(epsilon_bar_1){diag(x = 1 - epsilon_bar_1, nrow = 4) + 
  epsilon_bar_1 / 3 * (matrix(data = 1, nrow = 4, ncol = 4) - diag(x = 1, nrow = 4))}

#' @examples
#' sim_pi_igs(i=3,g = 3, s=5,alpha_pi=.1)

epsilon_iba_f <- function(i,epsilon_bar_1){ epsilon_ba=epsilon_ba_f(epsilon_bar_1)
epsilon_iba=plyr::raply(i,epsilon_ba)|>(`dimnames<-`)(list(i=1:i,b=1:4,a=1:4))
}


#' @params
#' @examples
#' sim_n_vsa(s = 3, n = 1000, v = 50, g = 5)
sim_n_vsa <- function(n=1000,# expeted sample size
                      tau_vgb = NULL,
                      pi_gs = NULL,
                      v = NULL,
                      g = NULL,
                      s = NULL,
                      epsilon_bar_1 = .001,
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
  epsilon_ba <- epsilon_ba_f(epsilon_bar_1)
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

#' @params
#' @examples
#' sim_n_vsa(s = 3, n = 1000, v = 50, g = 5)
sim_n_ivsa <- function(i,
                      n=1000,# expeted sample size
                      tau_ivgb = NULL,
                      pi_igs = NULL,
                      v = NULL,
                      g = NULL,
                      s = NULL,
                      epsilon_bar_1 = .001,
                      alpha_pi = 1) {


n_ivsa <- plyr::aaply(matrix(1:i)|>(`dimnames<-`)(list(i=1:i)),1,function(ii){
  sim_n_vsa(n = n, 
            tau_vgb = tau_ivgb[ii,,,,drop=FALSE]|>abind::adrop(1), 
            pi_gs = pi_igs[ii,,,drop=FALSE]|>abind::adrop(1),
            epsilon_bar_1 = epsilon_bar_1)})}


#' @examples
#' g = 5;v = 50; s = 3; n = 100; alpha_pi = .1
#' sim_tau_pi_epsilon_n(g = 5, v = 50, s = 3, n = 100, alpha_pi = 1)
sim_tau_pi_epsilon_n <- function(v, g, s, n, epsilon_bar_1 = .001, alpha_pi = 1) {
  tau_vgb <- sim_tau_vgb(v = v, g = g)
  pi_gs <- sim_pi_gs(g = g, s = s, alpha_pi = alpha_pi)
  reorder_g=order(pi_gs|>plyr::aaply(1,sum),decreasing = TRUE)|>
    (function(x){x|>setNames(x)})()
  pi_gs=pi_gs[reorder_g,]|>(`dimnames<-`)(list(g=1:g,s=1:s))
  tau_vgb=tau_vgb[,reorder_g,]|>(`dimnames<-`)(list(v=1:v,g=1:g,b=nucleotides))
  epsilon_ba=epsilon_ba_f(epsilon_bar_1)
  n_vsa <- sim_n_vsa(n = n, tau_vgb = tau_vgb, pi_gs = pi_gs, epsilon_bar_1 = epsilon_bar_1)
  list(tau_vgb = tau_vgb, pi_gs = pi_gs, epsilon_ba=epsilon_ba,n_vsa = n_vsa)
}





sim_tau_pi_epsilon_n_i <- function(i,v, g, s, n, epsilon_bar_1 = .001, alpha_pi = 1) {
  tau_ivgb <- sim_tau_ivgb(i=i,v = v, g = g)
  pi_igs <- sim_pi_igs(i=i,g = g, s = s, alpha_pi = alpha_pi)
  reorder_ig=plyr::aaply(pi_igs,1,function(x){order(x|>plyr::aaply(1,sum),decreasing = TRUE)|>
    (function(x){x|>setNames(x)})()})
  pi_igs=plyr::aaply(matrix(1:i)|>(`dimnames<-`)(list(i=1:i)),1,function(ii){pi_igs[ii,reorder_ig[ii,],]|>
                       (`dimnames<-`)(list(g=1:g,s=1:s))})
  tau_ivgb=plyr::aaply(matrix(1:i)|>(`dimnames<-`)(list(i=1:i)),
                       1,
                       function(ii){
                         tau_ivgb[ii,,reorder_ig[ii,],,drop=FALSE]})|>(`dimnames<-`)(list(i=1:i,v=1:v,g=1:g,b=nucleotides))
  epsilon_iba=epsilon_iba_f(i,epsilon_bar_1)
  
  rho_ivga=einsum::einsum('ivgb,iba->ivga',tau_ivgb,epsilon_iba)
  chi_ivsag=einsum::einsum('ivga,igs,iba->ivsag',rho_ivga,pi_gs)
  chi_ivsabg=einsum::einsum('ivgb,igs,iba->ivsabg',tau_ivgb,pi_igs,epsilon_iba)
  
  n_ivsa <- sim_n_ivsa(i,
                       n=n,# expeted sample size
                       tau_ivgb = tau_ivgb,
                       rho_ivga=rho_ivga,
                       chi_ivsag=chi_ivsag,
                       chi_ivsabg=chi_ivsabg,
                       pi_igs = pi_igs,
                       v = v,
                       g = g,
                       s = s,
                       epsilon_bar_1 = .001,
                       alpha_pi = 1)
  list(tau_ivgb = tau_vgb, pi_igs = pi_igs, epsilon_iba=epsilon_iba,n_vsa = n_ivsa,rho_ivga=rho_ivga,chi_ivsag=chi_ivsag,chi_ivsabg=chi_ivsabg)
}
