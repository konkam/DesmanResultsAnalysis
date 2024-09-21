reduce_tau_vgb_f<-function(tau_vgb,n_vsa,discriminant_v=discriminant_v_f(tau_vgb)){

X=tau_vgb[discriminant_v,,]|>
  aperm(c(1,3,2))|>
  array(dim=c(dim(tau_vgb[discriminant_v,,])[1]*
                dim(tau_vgb[discriminant_v,,])[3],
              dim(tau_vgb[discriminant_v,,])[2]))
Y=n_vsa|>
  aperm(c(1,3,2))|>
  array(dim=c(dim(n_vsa[discriminant_v,,])[1]*dim(n_vsa[discriminant_v,,])[3],dim(n_vsa)[2]))                 

plyr::alply(Y,2,function(y){nnls::nnls(X, y)})->
  nnls_outputs

nnls_outputs|>
  plyr::laply(`[[`,"x")|>
  (`>`)(0)|>
  plyr::aaply(2, any)|>
  which()->selected_g

nnls_outputs|>
  plyr::laply(`[[`,"x")|>
  plyr::aaply(1,`[`,selected_g)|>
  t()|>
  plyr::aaply(2,function(x){x/mean(x)})->pi_gs_0



list(tau_vgb=tau_vgb[,selected_g,,drop=FALSE],
     pi_gs_0=pi_gs_0)
}
