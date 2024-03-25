

if(FALSE){
  
  real_values_pi <- expand_grid(variant_index = 1:G, sample = 1:S) |>
    mutate(real_value = mapply(FUN = function(g, s) {
      pi_gs[g, s]
    }, variant_index, sample, SIMPLIFY = F) |> unlist())
  
  plot_pi_summary <- function(pi_summary, real_values_pi) {
    pi_summary |>
      mutate(variant_index = parname |>
               lapply(FUN = function(nm) {
                 nm |>
                   str_split_i(pattern = "\\[", i = 2) |>
                   str_split_i(pattern = ",", i = 1)
               }) |>
               unlist() |>
               as.numeric()) |>
      mutate(sample = parname |>
               lapply(FUN = function(nm) {
                 nm %>%
                   str_split_i(pattern = ",", i = 2) %>%
                   gsub("\\]", "", .)
               }) |>
               unlist() |>
               as.numeric()) |>
      left_join(real_values_pi) |>
      ggplot(aes(x = variant_index)) +
      theme_bw() +
      facet_wrap(~sample, ncol = 1, labeller = label_both) +
      geom_segment(aes(xend = variant_index, y = Lower95, yend = Upper95)) +
      geom_point(aes(y = Median), colour = "black") +
      geom_point(aes(y = real_value), colour = "red") +
      ggtitle("Real value in red, estimation in black") +
      scale_x_continuous(breaks = function(x) unique(x))
    # scale_x_discrete()
  }
  
  jags_samples_fixed_epsilon |>
    jags_sample_to_summary_tibble() |>
    plot_pi_summary(real_values_pi)
  
  
  
  #### Model code
  plot_jags_samples
  
  jags_samples |>
    jags_sample_to_summary_tibble() |>
    dplyr::filter(grepl("pi_gs", parname)) |>
    plot_pi_summary(real_values_pi)
  
  real_values_epsilon <- data.frame(a = 1:2) |>
    mutate(real_value = mapply(FUN = function(a) {
      tildeepsilon[a]
    }, a, SIMPLIFY = F) |> unlist()) |>
    mutate(parname = mapply(FUN = function(a) {
      paste("tildeepsilon[", a, "]", sep = "")
    }, a, SIMPLIFY = F) |> unlist())
  
  jags_samples$mcmc |>
    plyr::llply(`[`,,"tildeepsilon[2]")|>
    do.call(what=cbind)|>
    reshape2::melt()|>
    ggplot(aes(x = Var1,group=as.factor(Var2),col=as.factor(Var2),y=value)) +
    theme_bw() +
    geom_line()+
    scale_y_continuous(trans="log10")+
    geom_hline(yintercept = error_rate, colour = "red", size = 0.75) +
    facet_grid(~Var2)+
    ggtitle("Real value in red, estimation in black")
  
  jags_samples$mcmc |>
    plyr::llply(`[`,,"tildeepsilon[2]")|>
    do.call(what=cbind)|>
    reshape2::melt()|>
    ggplot(aes(x=as.factor(Var2),y=value)) +
    theme_bw() +
    geom_violin()+
    geom_point(size=.01,alpha=.5,position = position_jitter(w = 0.1, h = 0))+
    scale_y_continuous(trans="log10")+
    geom_hline(yintercept = error_rate, colour = "red", size = 0.75) +
    ggtitle("Real value in red, estimation in black")
  
  
}