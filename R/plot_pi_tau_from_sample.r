importance_g_it_f<-function(smc_samples,n_vsa,traces=NULL,
                          reorder_g=smc_reorder(smc_samples,n_vsa,traces=traces)){
  X=if(!is.null(smc_samples)){
    smc_samples|>
      mcmc_output_df(variable_name = "pi",index="gs")
  }else{convert_custom_ouput(traces["pi_tigs"])$pi_gs}
  X|>
    dplyr::group_by(iteration,chain,g)|>
    dplyr::summarise(part_of_g=sum(value))|>
    dplyr::ungroup()|>
    dplyr::left_join(reorder_g)
  }

convert_columns <- function(df) {
  df[] <- lapply(df, function(x) {
    if (is.factor(x)) {
      levels(x)[x]}else{x}})
    df[] <- lapply(df, function(x) {
      if (is.character(x)) {
        if (all(grepl("^[-+]?[0-9]*\\.?[0-9]+$", x))) {
        as.numeric(x)
      } else {
        x
      }
    } else {
      x
    }
  })
  return(df)
}


convert_custom_ouput<-function(traces){
  output=list()
  for(nm in setdiff(names(traces),"ess_t")){
    index=stringr::str_extract(nm, "[^_]+$")
    nom=stringr::str_extract(nm, ".*(?=_[^_]*$)")
    nom2=nom|>
      paste0("_", stringr::str_replace_all(index, "[ti]", ""))
    output[[nom2]]<-
      plyr::adply(
        traces[[nm]]|>
          as.array()|>
          namedims(index=index),1:length(dim(traces[[nm]]|>as.array())))|>
      convert_columns()|>
      dplyr::rename(chain=i,iteration=t,value="V1")|>
      dplyr::mutate(variable={nom})
    }
  output}



plot_pi_tau_from_sample<-
  function(smc_samples,
           traces=NULL,
           what=c("tau","pi","epsilon"),
           n_vsa,
           true_parameter=NULL,
           sel_v=1:2,
           sel_i=1:3,
           t_step=NULL,
           reorder_g=smc_reorder(smc_samples,n_vsa,traces=traces),
           importance_g_it=
             importance_g_it_f(smc_samples=smc_samples,
                               n_vsa=n_vsa,
                               reorder_g=reorder_g,traces=traces)){
    if(!is.null(true_parameter)){
      pi_g_true=true_parameter$pi_gs|>
        plyr::aaply(1,sum,.drop=FALSE)|>
        plyr::adply(1,.drop=FALSE)|>
        convert_columns()|>
        dplyr::rename(pi_g="V1")
      pi_gs_true=true_parameter$pi_gs|>
        plyr::adply(1:2,.drop=FALSE)|>
        convert_columns()|>
        dplyr::rename(value="V1")|>
        dplyr::group_by(g)|>
        dplyr::mutate(reorder_g=g,variable="pi",chain="true")|>
        dplyr::ungroup()|>
        (function(x){
          rbind(x|>
                  dplyr::mutate(iteration=0),
                x|>
                  dplyr::mutate(iteration=t_step))})()
     tau_true=true_parameter$tau_vgb|>plyr::adply(1:3)|>
       convert_columns()|>
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
    
    if(!is.null(traces)){traces_dfs<-convert_custom_ouput(traces)}
    if(is.element("tau",what)){
      
  (if(!is.null(smc_samples)){smc_samples|>
        mcmc_output_df(variable_name = "tau",index="vgb")|>
      dplyr::mutate(b=nucleotides[strtoi(b)])
    }else{traces_dfs$tau_vgb})|>
        dplyr::left_join(importance_g_it)->X
      if(!is.null(true_parameter)){
        X<-X|>
        (function(x){rbind(x,tau_true[names(x)])})()}
        
        X|>
          dplyr::filter(is.element(v,{{sel_v}}),
                        is.element(chain,{{sel_i}}),
                        is.null(t_step)|(strtoi(iteration)%%t_step==0))|>
          dplyr::mutate(reorder_g=paste0("g=",reorder_g),
                        v        =paste0("v=",v),
                        chain        =paste0("i=",chain))->X
        
        plot_tau=X|>
          #dplyr::mutate(g=paste0("g=",g),v=paste0("v=",v),chain=paste0("i=",i))|>
          ggplot(aes(x=iteration,y=value*part_of_g,group=interaction(chain,reorder_g,b,v),fill=b))+
          geom_bar(position = "stack", stat="identity",width=t_step)+
          xlab("")+ylab("")+
          theme(legend.position="none")+
          facet_grid(reorder_g+v~chain,scales="free_x")
        
        
        
        
        }else{plot_tau=NULL}  
      
      if(is.element("epsilon",what)){
      epsilon_tiba=if(!is.null(smc_samples)){
        if(any(grepl(pattern = "epsilon",smc_samples$mcmc[[1]]|>colnames()))){
          smc_samples|>
            mcmc_output_df(variable_name = "bar_epsilon[1]")
        }else{
            NULL}
      }else{
        if(!is.null(traces_dfs$epsilon_ba)){
          traces_dfs$epsilon_ba|>
            dplyr::filter(a=="a",b=="c")|>
            dplyr::mutate(value=3*value)}else{NULL}}
      
      plot_epsilon=if(is.null(epsilon_tiba)){
        epsilon_tiba|>
        ggplot(aes(x=iteration,y=value))+
        geom_line()+
        scale_y_continuous(trans="log10")+
        facet_wrap(~chain)}else{NULL}
      }else{plot_epsilon=NULL}
      
      if(is.element("pi",what)){
      pi_tigs=if(!is.null(smc_samples)){
        smc_samples|>
          mcmc_output_df(variable_name = "pi",index="gs")
      }else{
        if(!is.null(traces_dfs$pi_gs)){
          traces_dfs$pi_gs}else{NULL}}
      
      X=pi_tigs|>
       dplyr::left_join(reorder_g)
      
      if(!is.null(true_parameter)){X=X|>
       (`[`)(names(pi_gs_true))|>
       rbind(pi_gs_true)}
     plot_pi=   X|>ggplot(aes(x=iteration,y=value,group=reorder_g,fill=as.factor(reorder_g)))+
       geom_area()+
        facet_grid(s~chain,scales="free_x")
      }else{plot_pi=NULL}
    
      
      
      list(plot_tau=plot_tau,
           plot_epsilon=plot_epsilon,
           plot_pi=plot_pi)
}
