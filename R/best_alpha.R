#'@examples
#'g=200
#'alpha=.1
#'n=30000
#'selected_g_count(.1,g,n)
#'x=plyr::aaply(sort(outer(0:4,c(1,2,5),function(x,y){y*10^(-x)})),1,
#'selected_g_count,g=300,n=30000,.progress="text")
#'x


selected_g_count<-function(alpha,g,n){
  nrep=10000
  pis=LaplacesDemon::rdirichlet(nrep,
                            alpha = rep(alpha,g))
  pis_summary<-
    pis|>plyr::aaply(.margins = 1,.fun = function(p){
      c(seuil=1/n,
        npetit=sum(p<1/n),
        medianpetit=median(p[p<1/n]),
        ngrand=sum(p>1/n),
        mediangrand=median(p[p>1/n]))
    })|>
    plyr::aaply(2,median)

  c(pis_summary,
  plyr::aaply(pis,1,
                   function(p){sum(rmultinom(1,n,p)>0)})|>
    (function(x){c(alpha=alpha,min=min(x),median=median(x),mean=mean(x),max=max(x))})())}

selected_g_count_mean<-function(alpha,g,n){
  nrep=10000
  pis=LaplacesDemon::rdirichlet(nrep,
                                alpha = rep(alpha,g))
  pis<-pis[!is.na(plyr::aaply(pis,1,sum)),]
  
  plyr::aaply(pis,1,
              function(p){
                sum(rmultinom(1,n,p)>0)
              })|>
    mean(na.rm=TRUE)}



#'@examples
#'g=6
#'alpha=.1
#'n=30000
#'L(.1,g,5,n)
L<-function(alpha,g,g_target,n){
  (selected_g_count_mean(alpha,g,n)-g_target)^2
}

#'@examples
#'g=6
#'g_target=4
#'alpha=.1
#'n=30000
#'alpha=best_alpha(g,g_target,n)
best_alpha=function(g,g_target,n){
  alpha=sort(outer(0:4,c(1:9),function(x,y){y*10^(-x)}))
  ls=plyr::aaply(alpha,1,L,g=g,g_target=g_target,n=n,.progress="text")
  #'@examples
  #'library(ggplot2);ggplot(data.frame(x=alpha,y=ls),aes(x=x,y=y))+geom_line()+scale_x_continuous(trans="log10")
  i=which.min(ls)
  alpha[i]
  #optimize(f = L,interval = c(1e-4,10),g=g,g_target=g_target,n=n)
}

