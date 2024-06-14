smc_reorder<-function(smc_samples,n_vsa){
  n_vsa|>plyr::aaply(2,sum,.drop=FALSE)|>
    plyr::adply(1)|>
    dplyr::rename("n_s"=V1)->n_s
  
  smc_samples|>
    mcmc_output_df(variable_name = "pi",index="gs")|>
    dplyr::left_join(n_s)|>
    dplyr::mutate(n_gsit=value*n_s)|>
    dplyr::group_by(iteration,chain,g)|>
    dplyr::summarise(part_of_g=sum(n_gsit))|>
    dplyr::ungroup()|>
    dplyr::group_by(iteration,chain)|>
    dplyr::mutate(reorder_g=rank(-part_of_g,ties.method = "first"))|>
    dplyr::ungroup()|>dplyr::select(-part_of_g)
}
