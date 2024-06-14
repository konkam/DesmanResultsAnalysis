

jags_model_string_fixed_variants_f_old <- 
  function(bar_epsilon = NA) {
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
        mixed_variants[v, g, a] = inprod(tau_vgb[v,g,], epsilon[,a])
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
      if (is.na(bar_epsilon_1)) {
        "bar_epsilon~ddirch(shape_epsilon)
"
      },
      "
     for (a in 1:4){
      for (b in 1:4){
        epsilon[a,b] = (a!=b)*bar_epsilon[2]/3 +(a==b)*bar_epsilon[1]
        }
    }
}
"
    )}
