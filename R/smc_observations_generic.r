#### Prior specification for the error matrix
#' @description
#' Let's call \eqn{\tilde\epsilon} the error rate, which is of order 10e-3.
#' \deqn{\epsilon=(1-\tilde\epsilon)\times I_4+ \frac{\tilde\epsilon}3\times(J_4-I_4),  } where \eqn{(\tilde\epsilon\sim\mathrm{Beta}(a,b)} means that the prior expected value for the diagonal element \eqn{\tilde\epsilon=\epsilon_{1,1}} is \eqn{\frac{a}{a+b}} and its variance is \eqn{\frac{ab}{(a+b)^2(a+b+1)}}
#' Moment-matching specification for \eqn{a,b}:
#' \deqn{
#' \begin{align}
#' &\frac{a}{a+b} = 1 - E[\tilde\epsilon]\implies b = a\frac{E[\tilde\epsilon]}{1-E[\tilde\epsilon]} \\
#' &\text{Var}(\tilde\epsilon) = \frac{ab}{(a+b)^2(a+b+1)} = \frac{(1-\eta)\eta}{\frac{a}{1-\eta}+1} = \frac{(1-\eta)^2\eta}{a+1-\eta}\implies a = \eta - 1 + \frac{(1-\eta)^2\eta}{\text{Var}(\epsilon_{1,1})}
#' \end{align}
#' }
#' @examples
#' alpha_bar_epsilon_specification(.5,.25)
alpha_bar_epsilon_specification <- function(bar_epsilon_1_mean, bar_epsilon_1_std) {
  a <- bar_epsilon_1_mean - 1 + (1 - bar_epsilon_1_mean)^2 * bar_epsilon_1_mean / bar_epsilon_1_std^2
  c(a,a * bar_epsilon_1_mean / (1 - bar_epsilon_1_mean))
}




smc_observation_f <- 
  function(n_vsa,
           gs="jags",
           G,
           alpha_pi=NULL,
           alpha_epsilon=NULL,
           bar_epsilon_1_std=NULL,
           bar_epsilon_1_mean=NULL,
           alpha_bar_epsilon=NULL,
           bar_epsilon_1=NULL,
           tau_vgb=NULL,
           alpha_tau=NULL,
           kappa_rho=NULL,
           alpha_rho=NULL) {
    
    
    dimnames(n_vsa) <- lapply(dim(n_vsa), seq_len)
    
    observations=list(n_vsa=n_vsa,
           n_vs=n_vsa |> apply(MARGIN = c(1, 2), FUN = sum))
    if(!is.null(tau_vgb)){observations=c(observations,list(tau_vgb=tau_vgb))}
    if(!is.null(alpha_epsilon)){observations=c(observations,list(alpha_epsilon=alpha_epsilon,rep_alpha_epsilon=rep(alpha_epsilon,4)))}
    if(!is.null(alpha_bar_epsilon)){observations=c(observations,list(alpha_bar_epsilon=alpha_bar_epsilon))}else{
      if(!is.null(bar_epsilon_1_std)&!is.null(bar_epsilon_1_mean)){
        alpha_bar_epsilon=alpha_bar_epsilon_specification(
          bar_epsilon_1_std =bar_epsilon_1_std,
          bar_epsilon_1_mean=bar_epsilon_1_mean)
        observations=c(observations,list(alpha_bar_epsilon=alpha_bar_epsilon))}}
    if(!is.null(kappa_rho)){observations=c(observations,list(kappa_rho=kappa_rho))}
    if(!is.null(alpha_rho)){observations=c(observations,list(alpha_rho=alpha_rho,rep_alpha_rho=rep(alpha_rho,4)))}

    if(!is.null(alpha_pi)){observations=c(observations,
                                   list(alpha_pi=alpha_pi,
                                        rep_alpha_pi=rep(alpha_pi,G)))}
    observations
    
    
    constants=list(
      V = dim(n_vsa)[1],
      S = dim(n_vsa)[2],
      G=G,
      G4=4^G)
    if(!is.null(alpha_tau)){observations=c(observations,list(alpha_tau=alpha_tau))}
    
    if(!is.null(bar_epsilon_1)){constants=c(constants,
                                            list(bar_epsilon_1=bar_epsilon_1,bar_epsilon=c(bar_epsilon_1,1-bar_epsilon_1)))}
    
    
    if(gs!="nimble"){c(constants,observations)}else{list(constants=constants,observations=observations)}
    }
