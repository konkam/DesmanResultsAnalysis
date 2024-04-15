
#'@examples
#'tar_load(mcmc_output)
#'plot_tilde_epsilon(mcmc_output)
plot_tilde_epsilon<-function(mcmc_output,error_rate){
  mcmc_output|>
    mcmc_output_df(variable_name = "tildeepsilon[2]")|>
  ggplot(aes(x = iteration,group=as.factor(variable),col=as.factor(variable),y=value)) +
  theme_bw() +
  geom_line()+
  scale_y_continuous(trans="log10")+
  geom_hline(yintercept = error_rate, colour = "red", size = 0.75) +
  facet_grid(~variable)+
  ggtitle("Real value in red, estimation in black")}


#'@examples
#'tar_load(mcmc_output)
#'tar_load(ln_vsa)
#'n_vsa=ln_vsa[[1]]
#'variants=1:4
#'plot_pi(mcmc_output,n_vsa,variants)

plot_pi<-function(mcmc_output,n_vsa,variants){
  if(is.null(names(variants))){names(variants)<-1:length(variants)}
  names(variants)<-ifelse(is.na(names(variants)),1:length(variants),names(variants))
mcmc_output |>
    mcmc_output_df(variable_name = "pi",index = "gs")|>
    dplyr::mutate(s=if(!is.null(dimnames(n_vsa)[[2]])){dimnames(n_vsa)[[2]][strtoi(s)]}else{s})|>
    dplyr::mutate(g=names(variants)[strtoi(g)])|>
    dplyr::filter(iteration%%100==0)|>
    dplyr::group_by(s,g)|>
    dplyr::mutate(meanpi=mean(value))|>
    dplyr::ungroup()|>
    dplyr::group_by(s,iteration,chain)|>
    dplyr::mutate(valuemax=max(value[meanpi==max(meanpi)]))|>
    dplyr::ungroup()|>
    dplyr::group_by(s)|>
    dplyr::mutate(rankpi =rank(valuemax,ties.method="average"),
                  rankpi=rankpi/max(rankpi))|>
    dplyr::ungroup()|>
    dplyr::arrange(rankpi)|>
    ggplot2::ggplot(ggplot2::aes(x = rankpi,
               group=g,fill=g,y=value)) +
    ggplot2::theme_bw() +
    ggplot2::geom_area()+
    ggplot2::facet_wrap(~s)+
    ggplot2::coord_flip()
}

