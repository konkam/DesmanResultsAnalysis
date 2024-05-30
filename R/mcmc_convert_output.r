
#'@examples
#'variable_name="pi"
#'index="gs"
#' tau_pi_n <- sim_tau_pi_epsilon_n(v = 50, g = 5, s = 3, n = 1000, alpha_pi = 1)
#' block_tau=FALSE
#' gs="jags"
#' block_tau=FALSE
#' mcmc_output=mcmc_desman_run(burnin = 10,sample = 50,n.chains = 2,gs="jags",   n_vsa = tau_pi_n$n_vsa,   G = 5,block_tau=FALSE)
#'mcmc_output_df(mcmc_output,variable_name=NULL,index=NULL)
#'mcmc_output_df(mcmc_output,variable_name="tau",index="vsa)
mcmc_output_df<-function(mcmc_output,variable_name=NULL,index=NULL){
  mcmc_output$mcmc |>
    (function(x){x|>setNames(1:length(x))})()|>
    c(list(along=3))|>
    do.call(what=abind::abind)|>
    reshape2::melt()|>
    setNames(c("iteration","variable","chain","value"))->output_df
  if(!is.null(variable_name)){
    output_df<-output_df|>dplyr::filter(grepl(x = variable,pattern=variable_name,fixed = TRUE))
    
  }
  if(!is.null(index)){
    for(i in 1:nchar(index)){
      output_df[[substr(index,i,i)]]<-
        stringr::str_split_i(stringr::str_extract(output_df$variable,"[,0-9]+"),",",i)
    }
  }
  output_df
  }
#'@examples
#'#'variable_name="pi"
#'index="gs"
#' tau_pi_n <- sim_tau_pi_epsilon_n(v = 50, g = 5, s = 3, n = 1000, alpha_pi = 1)
#' block_tau=FALSE
#' gs="jags"
#' block_tau=FALSE
#' mcmc_output=mcmc_desman_run(burnin = 10,sample = 50,n.chains = 2,gs="jags",   n_vsa = tau_pi_n$n_vsa,   G = 5,block_tau=FALSE)
#'mcmc_output_array(mcmc_output,variable_name="pi",index="gs")
#'mcmc_output_array(mcmc_output,variable_name="tau",index="vsa")[]

mcmc_output_array<-function(mcmc_output,variable_name=NULL,index=NULL){
  mcmc_output$mcmc |>
    (function(x){x|>setNames(1:length(x))})()|>
    c(list(along=3))|>
    do.call(what=abind::abind)->output_array
  if(!is.null(variable_name)){
    output_array<-output_array[,grepl(variable_name,sort(dimnames(output_array)[[2]])),]
  }
  if(!is.null(index)){
    dimnames_x=list()
    dim_x=c()
    for(i in 1:nchar(index)){
      dimnames_x[[substr(index,i,i)]]<-
        stringr::str_split_i(stringr::str_extract(dimnames(output_array)[[2]],"[,0-9]+"),",",i)|>unique()
      dim_x[substr(index,i,i)]<-n_distinct(dimnames_x[[substr(index,i,i)]])
    }
    output_array=array(output_array,dim=c(dim(output_array)[1],dim_x,dim(output_array)[3]))|>
      (`dimnames<-`)(c(dimnames(output_array)[1],dimnames_x,dimnames(output_array)[3]))
  }
  output_array
}
