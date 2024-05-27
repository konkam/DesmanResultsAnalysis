
#'@examples
#'variable_name="pi"
#'index="gs"
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
  