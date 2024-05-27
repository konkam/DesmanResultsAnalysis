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
"R" |>
  list.files(full.names = TRUE) |>
  sapply(FUN = source)

list(
  tar_target(
    desman_input_file,
    "/work_projet/ala/metachick-fugace/Analyses_Ouléye/DESMAN_2023/DESMAN_iterAddPost/DESMAN_blaTEM_1A_1_HM749966/freq2Desman.txt"
  ),
  tar_target(
    fasta_file,
    "/work_projet/ala/metachick-fugace/BD/allTEM.fasta"
  ),
  tar_target(
    desman_tuning_parameters,
    list(
      sample = 10000, burnin = 1000, adapt = 1000,
      error_rate = 0.001,
      prior_std = 0.01,
      n_chains = 2,
      alpha0 = .1
    )
  ),
  tar_target(fixed_tildeepsilon, .0001),
  tar_target(
    name =
      variants,
    command =
      fasta_file |>
        get_data_from_server() |>
        read_fasta_file() |>
        (function(x) {
          x[nchar(x) == 861]
        })()
  ),
  tar_target(
    name = ln_vsa,
    command = desman_input_file |> get_data_from_server() |> read_desman_input_files()
  ),
  tar_target(
    name = tau_vga,
    command = variants |>
      translate_dna_string_vector_to_string_matrix() |>
      translate_dna_matrix_to_binary_array()
  ),
  tar_target(
    name = mcmc_output,
    command = do.call(
      what = mcmc_fixed_variants,
      args = c(
        list(
          n_vsa = ln_vsa[[1]],
          tau_vga = tau_vga,
          G = dim(tau_vga)[2] - 1,
          tildeepsilon = NA
        ),
        desman_tuning_parameters
      )
    )
  ),
  tar_target(
    name = mcmc_output_fixed_epsilon,
    command = do.call(
      what = mcmc_fixed_variants,
      args = c(
        list(
          n_vsa = ln_vsa[[1]],
          tau_vga = tau_vga,
          G = dim(tau_vga)[2] - 1,
          tildeepsilon = fixed_tildeepsilon
        ),
        desman_tuning_parameters
      )
    )
  ),
  tarchetypes::tar_render(
    name = report,
    path = "reports/fixed_variants.rmd",
    output_file = "fixed_variants_large_dictionnary.md"
  )
)
