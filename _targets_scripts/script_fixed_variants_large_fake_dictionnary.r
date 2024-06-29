# Load packages required to define the pipeline:
library(targets)
library(tidyverse)
library(ggplot2)
library(dplyr)
library(runjags)

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
"R"|>list.files(full.names = TRUE)|>sapply(FUN = source)->noprint
desman_input_file<-"/work_projet/ala/metachick-fugace/Analyses_Ouléye/DESMAN_2023/DESMAN_iterAddPost/DESMAN_erm_F__3_M17808/freq2Desman.txt"


list(
  tar_target(name=
               variants0,
             command=
               c("atgacaaaaaagaaattgcccgttcgttttacgggtcagcactttactattgataaagtgctaataaaagatgcaataagacaagcaaatataagtaatcaggatacggttttagatattggggcaggcaaggggtttcttactgttcatttattaaaaatcgccaacaatgttgttgctattgaaaacgacacagctttggttgaacatttacgaaaattattttctgatgcccgaaatgttcaagttgtcggttgtgattttaggaattttgcagttccgaaatttcctttcaaagtggtgtcaaatattccttatggcattacttccgatattttcaaaatcctgatgtttgagagtcttgGaaattttctgggagCttccattGtccttcaattagaacctacacaaaagttattttcgaggaagctttacaatccatataccgttttctatcatactttttttgatttgaaacttgtctatgaggtaggtcctgaaagtttcttgccaccgccaactgtcaaatcagccctgttaaacattaaaagaaaacacttattttttgattttaagtttaaagccaaatacttagcatttatttcctgtctgttagagaaacctgatttatctgtaaaaacagctttaaagtcgattttcaggaaaagtcaggtcaggtcaatttcggaaaaattcggtttaaaccttaatgctcaaattgtttgtttgtctccaagtcaatggttaaactgttttttggaaatgctggaagttgtccctgaaaaatttcatccttcgtag",
                 "atgacaaaaaagaaattgcccgttcgttttacgggtcagcactttactattgataaagtgctaataaaagatgcaataagacaagcaaatataagtaatcaggatacggttttagatattggggcaggcaaggggtttcttactgttcatttattaaaaatcgccaacaatgttgttgctattgaaaacgacacagctttggttgaacatttacgaaaattattttctgatgcccgaaatgttcaagttgtcggttgtgattttaggaattttgcagttccgaaatttcctttcaaagtggtgtcaaatattccttatggcattacttccgatattttcaaaatcctgatgtttgagagtcttgGaaattttctgggagGttccattGtccttcaattagaacctacacaaaagttattttcgaggaagctttacaatccatataccgttttctatcatactttttttgatttgaaacttgtctatgaggtaggtcctgaaagtttcttgccaccgccaactgtcaaatcagccctgttaaacattaaaagaaaacacttattttttgattttaagtttaaagccaaatacttagcatttatttcctgtctgttagagaaacctgatttatctgtaaaaacagctttaaagtcgattttcaggaaaagtcaggtcaggtcaatttcggaaaaattcggtttaaaccttaatgctcaaattgtttgtttgtctccaagtcaatggttaaactgttttttggaaatgctggaagttgtccctgaaaaatttcatccttcgtag",
                 "atgacaaaaaagaaattgcccgttcgttttacgggtcagcactttactattgataaagtgctaataaaagatgcaataagacaagcaaatataagtaatcaggatacggttttagatattggggcaggcaaggggtttcttactgttcatttattaaaaatcgccaacaatgttgttgctattgaaaacgacacagctttggttgaacatttacgaaaattattttctgatgcccgaaatgttcaagttgtcggttgtgattttaggaattttgcagttccgaaatttcctttcaaagtggtgtcaaatattccttatggcattacttccgatattttcaaaatcctgatgtttgagagtcttgAaaattttctgggagGttccattGtccttcaattagaacctacacaaaagttattttcgaggaagctttacaatccatataccgttttctatcatactttttttgatttgaaacttgtctatgaggtaggtcctgaaagtttcttgccaccgccaactgtcaaatcagccctgttaaacattaaaagaaaacacttattttttgattttaagtttaaagccaaatacttagcatttatttcctgtctgttagagaaacctgatttatctgtaaaaacagctttaaagtcgattttcaggaaaagtcaggtcaggtcaatttcggaaaaattcggtttaaaccttaatgctcaaattgtttgtttgtctccaagtcaatggttaaactgttttttggaaatgctggaagttgtccctgaaaaatttcatccttcgtag",
                 "atgacaaaaaagaaattgcccgttcgttttacgggtcagcactttactattgataaagtgctaataaaagatgcaataagacaagcaaatataagtaatcaggatacggttttagatattggggcaggcaaggggtttcttactgttcatttattaaaaatcgccaacaatgttgttgctattgaaaacgacacagctttggttgaacatttacgaaaattattttctgatgcccgaaatgttcaagttgtcggttgtgattttaggaattttgcagttccgaaatttcctttcaaagtggtgtcaaatattccttatggcattacttccgatattttcaaaatcctgatgtttgagagtcttgGaaattttctgggagGttccattAtccttcaattagaacctacacaaaagttattttcgaggaagctttacaatccatataccgttttctatcatactttttttgatttgaaacttgtctatgaggtaggtcctgaaagtttcttgccaccgccaactgtcaaatcagccctgttaaacattaaaagaaaacacttattttttgattttaagtttaaagccaaatacttagcatttatttcctgtctgttagagaaacctgatttatctgtaaaaacagctttaaagtcgattttcaggaaaagtcaggtcaggtcaatttcggaaaaattcggtttaaaccttaatgctcaaattgtttgtttgtctccaagtcaatggttaaactgttttttggaaatgctggaagttgtccctgaaaaatttcatccttcgtag")),
  tar_target(name=
               variants,
             command=
               c(variants0,create_fake_variants(n = 12,dic=variants0,taux = .05))),
  tar_target(name =ln_vsa,
    command = desman_input_file|>get_data_from_server()|>read_desman_input_files()),
  tar_target(
    name=tau_vgb,
    command =variants|>
      translate_dna_string_vector_to_string_matrix()|>
      translate_dna_matrix_to_binary_array()),
  tar_target(name = mcmc_output,
             command = mcmc_fixed_variants(ln_vsa[[1]][,1,,drop=FALSE],
                                                       tau_vgb,
                                                       G=dim(tau_vgb)[2]-1,
                                             bar_epsilon=NA,
                                                       bar_epsilon_1 = 0.001,
                                                       prior_std = 0.01,
                                                       n_chains=2,
                                                       alpha_pi=.1,
                                             burnin = 40,
                                             sample = 1000,
                                             adapt=500)),
  tar_target(name = mcmc_output_fixed_epsilon,
             command = mcmc_fixed_variants(ln_vsa[[1]][,1,,drop=FALSE],
                                             tau_vgb,
                                             G=dim(tau_vgb)[2]-1,
                                             bar_epsilon=.0001,
                                             bar_epsilon_1 = 0.001,
                                             prior_std = 0.01,
                                             n_chains=2,
                                             alpha_pi=.1)))
