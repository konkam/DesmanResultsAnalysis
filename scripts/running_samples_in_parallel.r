#load the rproject at the root of the Desman folder


# Run the R scripts in the R/ folder with your custom functions:
"R" |>list.files(full.names = TRUE) |>sapply(FUN = source)->noprint

E_coli_MG1655_mdh_C995_zip= "inst/extdata/script_20240709/E_coli_MG1655_mdh_C995.zip"
fastas=E_coli_MG1655_mdh_C995_zip|>
  utils::unzip(exdir=tempdir())|>
  sapply(read_fasta_file)
desman_input_file="inst/extdata/script_20240709/freq2Desman_E_coli_MG1655_mdh.txt"
fasta_file="inst/extdata/script_20240709/mdh_alleles_1size.fasta"
n_vsa= desman_input_file |> 
  get_data_from_server() |> 
  read_desman_input_files(s_pattern="[[:alnum:]]+")|>
  (`[[`)(1)
variants=
  fasta_file |>
  get_data_from_server() |>
  read_fasta_file()
tau_vgb=
  command = variants |>
  translate_dna_string_vector_to_string_matrix() |>
  translate_dna_matrix_to_binary_array()


g=5
t_max=4
n_chains=2


parallel_sample

library(parallel)
detectCores()
n_vsa|>
  plyr::alply(2,.drop=FALSE,identity())|>
  parallel::mclapply(FUN=function(n_vsa_one_sample,t_max,g,n_chains){
    smc_output_desman= desman_run(n_vsa=n_vsa_one_sample,
                                  t_max=t_max,
                                  g=g,
                                  n_chains =n_chains)
    
    smc_output_AM1= smc_custom(n_vsa=ln_vsa[[1]],
                               gs="custom",
                               g=g,
                               tau_vgb = NULL,
                               tau_vgb_0=NULL,
                               block_tau=FALSE,
                               bar_epsilon_1_std=NULL,
                               bar_epsilon_1_mean=NULL,
                               alpha_bar_epsilon=NULL,
                               bar_epsilon_1=.001,
                               alpha_rho=NULL,
                               alpha_pi=1,
                               n_chains = 3,
                               n_vsa_df=NULL,
                               alpha_tau=.001,
                               t_min=1,
                               t_max=4,
                               ess_min=NULL,
                               trace_all=TRUE,
                               mcmc=TRUE,
                               tempering_n=FALSE,
                               bar_epsilon_1_tempering=FALSE,
                               alpha_tau_tempering=FALSE,
                               alpha_tau_seq=NULL,
                               bar_epsilon_seq=NULL)
    #  tar_target(
    #    name = smc_output_AM1_block,
    #    command = smc_run_block_tau(
    #                n_vsa=ln_vsa[[1]],
    #                G=5,
    #                alpha_pi=.1,
    #                alpha_bar_epsilon=c(1,10),
    #                n_chains =3,
    #                mcmc=TRUE)
    #  ),
    smc_output_AM2_fixed_tau = smc_run_fixed_tau(n_vsa=ln_vsa[[1]],
                                                 tau_vgb = tau_vgb,
                                                 block_tau=FALSE,
                                                 alpha_pi=.1,
                                                 g=10,
                                                 alpha_bar_epsilon=c(1,10),
                                                 n_chains =n_chains,
                                                 mcmc=TRUE,
                                                 t_max=4)
    smc_output_AM3_relax_rho =  smc_run_relax_rho(n_vsa=ln_vsa[[1]],
                                                  g=5,
                                                  tau_vgb,
                                                  alpha_pi=.1,
                                                  alpha_rho=.1,
                                                  n_chains =n_chains,
                                                  mcmc=TRUE,
                                                  t_max=t_max) 
    
    smc_output_AM4_relax_tau=smc_run_relax_tau(n_vsa=ln_vsa[[1]],
                                               g=5,
                                               alpha_pi=.1,
                                               alpha_tau=.001,
                                               alpha_bar_epsilon=c(1,100),
                                               n_chains =n_chains,
                                               mcmc=t_max,
                                               t_max=4)
    
    list(smc_output_desman=smc_output_desman,
         smc_output_AM1=smc_output_AM1,you
         smc_output_AM2_fixed_tau=smc_output_AM2_fixed_tau,
         smc_output_AM3_relax_rho=smc_output_AM3_relax_rho,
         smc_output_AM4_relax_tau=smc_output_AM4_relax_tau)
    }
t_max=t_max,g=g,n_chains=n_chains)
