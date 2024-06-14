importance_g_it_f<-function(smc_samples,n_vsa,
                          reorder_g=smc_reorder(smc_samples,n_vsa)){
  smc_samples|>
    mcmc_output_df(variable_name = "pi",index="gs")|>
    dplyr::group_by(iteration,chain,g)|>
    dplyr::summarise(part_of_g=sum(value))|>
    dplyr::ungroup()|>
    dplyr::left_join(reorder_g)}


plot_pi_tau_from_sample<-
  function(smc_samples,
           n_vsa,
           true_parameter=NULL,
           v=1:2,
           t_step=NULL,
           reorder_g=smc_reorder(smc_samples,n_vsa),
           importance_g_it=
             importance_g_it_f(smc_samples=smc_samples,
                               n_vsa=n_vsa,
                               reorder_g=reorder_g)){

    if(!is.null(true_parameter)){
      pi_g_true=true_parameter$pi_gs|>
        plyr::aaply(1,sum,.drop=FALSE)|>
        plyr::adply(1,.drop=FALSE)|>
        dplyr::rename(pi_g="V1")
     tau_true=true_parameter$tau_vgb|>plyr::adply(1:3)|>
       dplyr::left_join(pi_g_true)|>
        dplyr::mutate(part_of_g=V1*pi_g,
                      reorder_g=g,
                      value=V1,
                      variable="tau",
                      chain="true")|>
       (function(x){
         rbind(x|>
           dplyr::mutate(iteration=0),
           x|>
           dplyr::mutate(iteration=t_step))})()
        
    }
      smc_samples|>
        mcmc_output_df(variable_name = "tau",index="vgb")|>
        dplyr::mutate(b=nucleotides[strtoi(b)])|>
        dplyr::left_join(importance_g_it)->X
      if(!is.null(true_parameter)){X<-X|>
        (function(x){rbind(x,tau_true[names(x)])})()}
      plot_tau=X|>
        dplyr::filter(is.element(v,{{v}}),
                      is.null(t_step)|(strtoi(iteration)%%t_step==0))|>
        ggplot(aes(x=iteration,y=value*part_of_g,group=interaction(chain,reorder_g,b,v),fill=b))+
        geom_area()+
        facet_grid(strtoi(reorder_g)+strtoi(v)~chain,scales="free_x")
      plot_epsilon=smc_samples|>
        mcmc_output_df(variable_name = "bar_epsilon[1]")|>
        ggplot(aes(x=iteration,y=value))+
        geom_line()+
        scale_y_continuous(trans="log10")+
        facet_wrap(~chain)
        
      
     plot_pi= smc_samples|>
        mcmc_output_df(variable_name = "pi",index="gs")|>
       dplyr::left_join(reorder_g)|>
        ggplot(aes(x=iteration,y=value,group=chain,colour=chain))+
        facet_grid(reorder_g~s)
      
      
      list(plot_tau=plot_tau,plot_epsilon=plot_epsilon,
           plot_pi=plot_pi)
}
