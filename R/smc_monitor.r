
smc_monitor_f <- 
  function(gs="jags",
           fixed_tau=FALSE,
           fixed_bar_epsilon=FALSE,
           constrained_epsilon_matrix=TRUE,
           relax_rho=FALSE,
           fixed_alpha_rho=FALSE){
    c("pi_gs",
      if(!fixed_tau&!relax_rho){"tau_vgb"},
      if(!fixed_bar_epsilon&constrained_epsilon_matrix){"bar_epsilon[1]"},
      if(!fixed_bar_epsilon&!constrained_epsilon_matrix){"epsilon_ba"},
      if(relax_rho){"rho_vga"},
      if(relax_rho&!fixed_alpha_rho){"alpha_rho"})}