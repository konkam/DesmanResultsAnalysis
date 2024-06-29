rdirichlet_smallalpha<-function (n = 1, alpha) {
  the_dim=length(alpha)
  Gam <- t(plyr::aaply(alpha,1,function(x){rgamma(n = n, shape = x)}))
  zeros<-plyr::aaply(Gam==0,1,all)
    if(any(zeros)){
    
    Gam[zeros,] <- t(rmultinom(sum(zeros),prob=rep(1/the_dim,the_dim),size =1))
    
  
  }
  Gam/rowSums(Gam)}

  
  
  
