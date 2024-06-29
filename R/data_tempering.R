
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)|>attach()
#'n_vsa_df<-reorder_reads(n_vsa)
#'

reorder_reads<-function(n_vsa,seed=1){
  set.seed(seed)
  n_vsa_df<-n_vsa|>
    plyr::adply(1:3,.fun=function(x){data.frame(V1=rep(1,x))})|>
    dplyr::mutate(rank=rank(runif(dplyr::n())))|>
    dplyr::arrange(rank)|>
    dplyr::mutate(lambda_threshold=dplyr::row_number()/dplyr::n())|>
    dplyr::select(-rank,-V1)}


#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)|>attach()
#'n_vsa_df<-reorder_reads(n_vsa)
#'n_vsa_lambda=sub_sample_counts(n_vsa_df,lambda=.5)
#'sub_sample_counts(n_vsa_df,lambda=0)
sub_sample_counts<-function(n_vsa_df,lambda){
  n_vsa_df|>
#    dplyr::filter(lambda_threshold<=lambda)|>
    #    dplyr::select(-lambda_threshold)|>
plyr::daply(.drop_i=FALSE,
                .drop_o=TRUE,
                .variables=~v+s+a,
                .fun=nrow)}
#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1;lambda=.5
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)|>attach()
#'n_vsa_df<-reorder_reads(n_vsa)
#'tempering_n_srs(n_vsa,lambda)

tempering_n_srs<-function(n_vsa,lambda,seed=1,n_vsa_df=reorder_reads(n_vsa,lambda)){
  sub_sample_counts(n_vsa_df,lambda)}


#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1;lambda=.09
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)|>attach()
#'n_vsa_df<-reorder_reads(n_vsa)
#'
#'https://www.census.gov/content/dam/Census/library/working-papers/2014/adrm/rrs2014-07.pdf
reorder_counts<-function(n_vsa,seed=1){
  set.seed(seed)
  H=prod(dim(n_vsa))
  data.frame(n_vsa_v=c(n_vsa),
             i=1:H,
             r=(1:H)[order(runif(H))])
}


#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1;lambda=.09
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)|>attach()
#'n_vsa_df<-reorder_counts(n_vsa,seed=1)
#'tempering_n_stratified(n_vsa,n_vsa_df=n_vsa_df,lambda=lambda)
tempering_n_stratified<-function(n_vsa,n_plus=sum(n_vsa),seed=1,n_vsa_df=reorder_counts(n_vsa,seed=seed),lambda){
  sample_size=ceiling(n_plus*lambda)
  df=n_vsa_df|>
    dplyr::mutate(n_vsa_lambda=n_vsa_v*lambda,
                  allocation0=floor(n_vsa_lambda),
                  left=(n_vsa_v-allocation0/lambda),
                  loss1=-2*lambda+(2*allocation0+1)/n_vsa_v,
                  loss2=1/n_vsa_v)|>
    dplyr::arrange(loss1,loss2,r)
  remaining=sample_size-sum(df$allocation0)
  H=nrow(n_vsa_df)
  while(remaining>0){
    df[1,"allocation0"]<-    df[1,"allocation0"]+1
    loss=df[1,"loss1"]      <-    if( df[1,"allocation0"]>=df[1,"n_vsa_v"]){Inf }else{df[1,"loss2"]}
    i=1
    while(df[i+1,"loss1"]<loss&i<H){
      i=i+1
    }
    #print(i)
    if(i>1){df[1:i,]=df[c(2:i,1),]}
    
    remaining=remaining-1  
    
  }
  
  df$allocation0[df$i]|>array(dim = dim(n_vsa),dimnames=dimnames(n_vsa))
}


#'@examples
#'n=1000;v=20;g=5;s=3;alpha_pi=.1
#'sim_tau_pi_epsilon_n(v=v, g=g, s=s, n=n, bar_epsilon_1 = .001, alpha_pi=alpha_pi)|>attach()
#'n_vsa_df<-reorder_reads(n_vsa)
#'unaccounted_for=unaccounted_for_f(n_vsa_df,lambda=.5)

unaccounted_for_f<-function(n_vsa_df,lambda){
  n_vsa_df|>dplyr::filter(lambda_threshold>lambda)}
