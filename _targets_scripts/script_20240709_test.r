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
    name = smc_output_,
    command = smc_custom(n_vsa=ln_vsa[[1]],
                         gs="custom",
                         G=5,
                         tau_vgb = NULL,
                         tau_vgb_0=NULL,
                         block_tau=FALSE,
                         bar_epsilon_1_std=.01,
                         bar_epsilon_1_mean=.001,
                         alpha_bar_epsilon=NULL,
                         bar_epsilon_1=NULL,
                         alpha_rho=NULL,
                         alpha_pi=1,
                         n_chains = 3,
                         n_vsa_df=NULL,
                         alpha_tau=.001,
                         t_min=1,
                         t_max=30,
                         ess_min=NULL,
                         trace_all=TRUE,
                         mcmc=TRUE,
                         tempering_n=FALSE,
                         bar_epsilon_1_tempering=FALSE,
                         alpha_tau_tempering=FALSE,
                         alpha_tau_seq=NULL,
                         bar_epsilon_seq=NULL)
  )
)

