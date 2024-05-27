tau_vga=c("aaaaaa","aaaccc","cccaaa","cccttt","cccggg")|>
  translate_dna_string_vector_to_string_matrix()|>
  translate_dna_matrix_to_binary_array()

tau_vga1=c("accaaa","caaggg","caaatt","caataa","accccc")|>
  translate_dna_string_vector_to_string_matrix()|>
  translate_dna_matrix_to_binary_array()

pi_gs<-matrix(c(.2,.3,.1,.1,.2),5,1)


n_vsa=sim_n_vsa(tau_vga=tau_vga,n = 10000,pi_gs = pi_gs,error_rate = .0001)
drop(n_vsa)
g=dim(tau_vga)[2]
s=1


plyr::raply(10,sampler_tau(tau_vga,pi_gs,epsilon_ba,n_vsa,
                      v=dim(tau_vga)[1],
                      g=dim(pi_gs)[1],
                      g_neq_g=g_neq_g_f(g))|>unname()|>
              translate_dna_binary_array_to_string_vector())|>rbind(
tau_vga|>
  translate_dna_binary_array_to_string_vector()|>unname())



plyr::raply(10,sampler_tau(tau_vga1,pi_gs,epsilon_ba,n_vsa,
                           v=dim(tau_vga)[1],
                           g=dim(pi_gs)[1],
                           g_neq_g=g_neq_g_f(g))|>unname()|>
              translate_dna_binary_array_to_string_vector())|>rbind(
                tau_vga1|>
                  translate_dna_binary_array_to_string_vector()|>unname())
  
  


#################################################################################

tau_vga=c("ac","ca")|>
  translate_dna_string_vector_to_string_matrix()|>
  translate_dna_matrix_to_binary_array()


rdirichlet(c(100,100),1)



#################################################################################

tau_vga=c("aa","cc","ac","ca")|>
  translate_dna_string_vector_to_string_matrix()|>
  translate_dna_matrix_to_binary_array()

pi_gs<-matrix(c(.4,.4,.1,.1),4,1)
pi_gs0<-matrix(c(.1,.1,.4,.4),4,1)


epsilon_ba <- epsilon_ba_f(error_rate)
alpha_g=.0001
n_vsa=sim_n_vsa(tau_vga=tau_vga,n = 10000,pi_gs = pi_gs,error_rate = .0001)

  xi_vsabg=sampler_xi(tau_vga=tau_vga,
                      pi_gs=pi_gs,
                      epsilon_ba = epsilon_ba)
  nu_vsab<-nu_from_xi(xi=xi_vsabg)
  mu_vsab<-mu_from_xi(xi=xi_vsabg)
  pi_gs1<-sampler_pi(mu_vsag = mu_vsab,alpha_g = alpha_g)
  
