display_n_vsa<-function(n_vsa){
  n_vsa|>aperm(c(1,3:2))|>
    array(c(dim(n_vsa)[1],prod(dim(n_vsa)[2:3])))|>
    `dimnames<-`(list(v=dimnames(n_vsa)[[1]],
                  `s-a`=paste0("s", t(outer(dimnames(n_vsa)[[2]],nucleotides,paste,sep="-")))))
}