
model.string.intractable.f <- function(tildeepsilon = NA, G, tau_vga) {
  Gd <- dim(tau_vga)[2]
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
    for (g in 1:Gd){
      for (a in 1:4){
        mixed_variants[v, g, a] = inprod(tau_vga[v,g,], epsilon[,a])
      }
    }
  }

  # Latent multinomial observation probability
  for (v in 1:V){
    for (s in 1:S){
      for (g in 1:Gd){
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
  for (s in 1:S){",
    if (G < Gd) {
      "
          pi_sel_gs[1:G, s] ~ ddirch(alpha[1:G])
          selected_[1:Gd, s] ~ dsample(rep(1,Gd),G)
          selected_pos=order(selected_[1:Gd, s]*(1:N))[1:G]
          non_selected_pos=order((rep(1,Gd)-selected_[1:Gd, s])*(1:N))[1:Gd-G]
          pi_gs[selected_pos, s]  =pi_sel_gs[1:G, s]
          pi_gs[non_selected_pos, s]  =rep(0,(Gd-G))
"
    } else {
      "
    pi_gs[1:Gd, s] ~ ddirch(alpha[1:Gd])
"
    },
    "}",
    if (is.na(tildeepsilon)) {
      "tildeepsilon~ddirch(c(aa, bb))"
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
}