# Load packages required to define the pipeline:
library(targets)
library(tarchetypes)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(runjags)
library(RColorBrewer)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tibble"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multicore")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
"R" |>list.files(full.names = TRUE) |>sapply(FUN = source)->noprint


list(
  tar_target(E_coli_MG1655_mdh_C995_zip, "inst/extdata/script_20240709/E_coli_MG1655_mdh_C995.zip",format = "file"),
  tar_target(fastas, E_coli_MG1655_mdh_C995_zip|>
               utils::unzip(exdir=tempdir())|>
              sapply(read_fasta_file)),
  tar_target(desman_input_file, "inst/extdata/script_20240709/freq2Desman_E_coli_MG1655_mdh.txt",format = "file"),
  tar_target(fasta_file, "inst/extdata/script_20240709/mdh_alleles_1size.fasta",format = "file"),
  tar_target(
    name = ln_vsa,
    command = desman_input_file |> 
      get_data_from_server() |> 
      read_desman_input_files(s_pattern="[[:alnum:]]+")),
  tar_target(
    name =
      variants,
    command =
      fasta_file |>
      get_data_from_server() |>
      read_fasta_file()),
  tar_target(
    name = tau_vgb,
    command = variants |>
      translate_dna_string_vector_to_string_matrix() |>
      translate_dna_matrix_to_binary_array()),
  tar_target(
    name = smc_output_desman,
    command = desman_run(n_vsa=ln_vsa[[1]][1:3,,],
                         t_max=4,
                                     g=5,
                                     n_chains =3)),
  tar_target(
    name = smc_output_AM1,
    command = smc_custom(n_vsa=ln_vsa[[1]][1:3,,],
                         gs="custom",
                         g=5,
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
  ),
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
  tar_target(
    name = smc_output_AM2_fixed_tau,#fixed tau
    command = smc_run_fixed_tau(n_vsa=ln_vsa[[1]],
                         g=1455,
                         tau_vgb = tau_vgb,
                         block_tau=FALSE,
                         alpha_pi=.1,
                         alpha_bar_epsilon=c(1,10),
                         n_chains =3,
                         mcmc=TRUE,
                         t_max=4)),
  tar_target(
    name = smc_output_AM3_relax_rho,
    command =  smc_run_relax_rho(n_vsa,
                           g=5,
                           tau_vgb,
                           alpha_pi=.1,
                           kappa_rho=c(1,100),
                           n_chains =3,
                           mcmc=TRUE,
                           t_max=4) 
  ),
  tar_target(
    name = smc_output_AM4_relax_tau,
    command = smc_run_relax_tau(n_vsa=ln_vsa[[1]][1:3,,],
                                g=5,
                                alpha_pi=.1,
                                alpha_tau=.001,
                                alpha_bar_epsilon=c(1,100),
                                n_chains =3,
                                mcmc=TRUE,
                                t_max=4)
  )
)

