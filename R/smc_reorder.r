smc_reorder<-function(smc_samples,n_vsa,traces=NULL){
  n_vsa|>plyr::aaply(2,sum,.drop=FALSE)|>
    plyr::adply(1)|>
    convert_columns()|>
    dplyr::rename("n_s"=V1)->n_s
  
  X=if(!is.null(smc_samples)){
    smc_samples|>
      mcmc_output_df(variable_name = "pi",index="gs")
  }else{convert_custom_ouput(traces["pi_tigs"])$pi_gs}
  X|>dplyr::left_join(n_s)|>
    dplyr::mutate(n_gsit=value*n_s)|>
    dplyr::group_by(iteration,chain,g)|>
    dplyr::summarise(part_of_g=sum(n_gsit))|>
    dplyr::ungroup()|>
    dplyr::group_by(iteration,chain)|>
    dplyr::mutate(reorder_g=rank(-part_of_g,ties.method = "first"))|>
    dplyr::ungroup()|>dplyr::select(-part_of_g)
}
