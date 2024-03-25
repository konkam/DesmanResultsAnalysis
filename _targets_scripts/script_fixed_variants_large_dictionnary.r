# Load packages required to define the pipeline:
library(targets)
library(ggplot2)
library(dplyr)
library(tidyverse)
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
"R"|>list.files(full.names = TRUE)|>sapply(FUN = source)
desman_input_file<-"/work_projet/ala/metachick-fugace/Analyses_Ouléye/DESMAN_2023/DESMAN_iterAddPost/DESMAN_erm_F__3_M17808/freq2Desman.txt"


list(
  tar_target(names=
               variants,
             command=
               c("atgacaaaaaagaaattgcccgttcgttttacgggtcagcactttactattgataaagtgctaataaaagatgcaataagacaagcaaatataagtaatcaggatacggttttagatattggggcaggcaaggggtttcttactgttcatttattaaaaatcgccaacaatgttgttgctattgaaaacgacacagctttggttgaacatttacgaaaattattttctgatgcccgaaatgttcaagttgtcggttgtgattttaggaattttgcagttccgaaatttcctttcaaagtggtgtcaaatattccttatggcattacttccgatattttcaaaatcctgatgtttgagagtcttgGaaattttctgggagCttccattGtccttcaattagaacctacacaaaagttattttcgaggaagctttacaatccatataccgttttctatcatactttttttgatttgaaacttgtctatgaggtaggtcctgaaagtttcttgccaccgccaactgtcaaatcagccctgttaaacattaaaagaaaacacttattttttgattttaagtttaaagccaaatacttagcatttatttcctgtctgttagagaaacctgatttatctgtaaaaacagctttaaagtcgattttcaggaaaagtcaggtcaggtcaatttcggaaaaattcggtttaaaccttaatgctcaaattgtttgtttgtctccaagtcaatggttaaactgttttttggaaatgctggaagttgtccctgaaaaatttcatccttcgtag",
                 "atgacaaaaaagaaattgcccgttcgttttacgggtcagcactttactattgataaagtgctaataaaagatgcaataagacaagcaaatataagtaatcaggatacggttttagatattggggcaggcaaggggtttcttactgttcatttattaaaaatcgccaacaatgttgttgctattgaaaacgacacagctttggttgaacatttacgaaaattattttctgatgcccgaaatgttcaagttgtcggttgtgattttaggaattttgcagttccgaaatttcctttcaaagtggtgtcaaatattccttatggcattacttccgatattttcaaaatcctgatgtttgagagtcttgGaaattttctgggagGttccattGtccttcaattagaacctacacaaaagttattttcgaggaagctttacaatccatataccgttttctatcatactttttttgatttgaaacttgtctatgaggtaggtcctgaaagtttcttgccaccgccaactgtcaaatcagccctgttaaacattaaaagaaaacacttattttttgattttaagtttaaagccaaatacttagcatttatttcctgtctgttagagaaacctgatttatctgtaaaaacagctttaaagtcgattttcaggaaaagtcaggtcaggtcaatttcggaaaaattcggtttaaaccttaatgctcaaattgtttgtttgtctccaagtcaatggttaaactgttttttggaaatgctggaagttgtccctgaaaaatttcatccttcgtag",
                 "atgacaaaaaagaaattgcccgttcgttttacgggtcagcactttactattgataaagtgctaataaaagatgcaataagacaagcaaatataagtaatcaggatacggttttagatattggggcaggcaaggggtttcttactgttcatttattaaaaatcgccaacaatgttgttgctattgaaaacgacacagctttggttgaacatttacgaaaattattttctgatgcccgaaatgttcaagttgtcggttgtgattttaggaattttgcagttccgaaatttcctttcaaagtggtgtcaaatattccttatggcattacttccgatattttcaaaatcctgatgtttgagagtcttgAaaattttctgggagGttccattGtccttcaattagaacctacacaaaagttattttcgaggaagctttacaatccatataccgttttctatcatactttttttgatttgaaacttgtctatgaggtaggtcctgaaagtttcttgccaccgccaactgtcaaatcagccctgttaaacattaaaagaaaacacttattttttgattttaagtttaaagccaaatacttagcatttatttcctgtctgttagagaaacctgatttatctgtaaaaacagctttaaagtcgattttcaggaaaagtcaggtcaggtcaatttcggaaaaattcggtttaaaccttaatgctcaaattgtttgtttgtctccaagtcaatggttaaactgttttttggaaatgctggaagttgtccctgaaaaatttcatccttcgtag",
                 "atgacaaaaaagaaattgcccgttcgttttacgggtcagcactttactattgataaagtgctaataaaagatgcaataagacaagcaaatataagtaatcaggatacggttttagatattggggcaggcaaggggtttcttactgttcatttattaaaaatcgccaacaatgttgttgctattgaaaacgacacagctttggttgaacatttacgaaaattattttctgatgcccgaaatgttcaagttgtcggttgtgattttaggaattttgcagttccgaaatttcctttcaaagtggtgtcaaatattccttatggcattacttccgatattttcaaaatcctgatgtttgagagtcttgGaaattttctgggagGttccattAtccttcaattagaacctacacaaaagttattttcgaggaagctttacaatccatataccgttttctatcatactttttttgatttgaaacttgtctatgaggtaggtcctgaaagtttcttgccaccgccaactgtcaaatcagccctgttaaacattaaaagaaaacacttattttttgattttaagtttaaagccaaatacttagcatttatttcctgtctgttagagaaacctgatttatctgtaaaaacagctttaaagtcgattttcaggaaaagtcaggtcaggtcaatttcggaaaaattcggtttaaaccttaatgctcaaattgtttgtttgtctccaagtcaatggttaaactgttttttggaaatgctggaagttgtccctgaaaaatttcatccttcgtag",
                 fake_1="atgacaaataagaaattgcccgttcgttttacgcgtcagcactttactattgataaagtgctaataaaagatgcaataagacaagcaaatataagtaatcaggatacggttttagatattggggcaggcaaggggtttcttactgttcatttattaaaaatcgccaacaatgttgttgctattgaaaacgacacagctttggttgaacatttacgaaaattattttctgatgcccgaaatgttcaagttgtcggttgtgattttaggaattttgcagttccgaaatttcctttcaaagtggtgtcaaatattccttatggcattacttccgatattttcaaaatcctgatgtttgagagtcttgGaaattttctgggagGttccattAtccttcaattagaacctacacaaaagttattttcgaggaagctttacaatccatataccgttttctatcatactttttttgatttgaaacttgtctatgaggtaggtcctgaaagtttcttgccaccgccaactgtcaaatcagccctgttaaacattaaaagaaaacacttattttttgattttaagtttaaagccaaatacttagcatttatttcctgtctgttagagaaacctgatttatctgtaaaaacagctttaaagtcgattttcaggaaaaatcaggtcaggtcaatttcggaaatattcggtttaaaccttaatgctcaaatagtttgtttgtctccaagtcaatggttaaactgtcttttggaaatgctggaagttgtgcctgaataatttcatccttcgtag",
                 fake_2="ttttagatattggggcaggcaaggggtttcttactgttcatttattaaaaatcgccaacaatgttgttgctattgaaaacgacacagctttggttgaacatttacgaaaattattttctgatgcccgaaatgttcaagttgtcggttgtgattttaggaattttgcagttccgaaatttcctttcaaagtggtgtcaaatattccttatggcattacttccgatattttcaaaatcctgatgtttgagagtcttgGaaattttctgggagGttccattAtccttcaattagaacctacacaaaagttattttcgaggaagctttacaatccatataccgttttctatcatactttttttgatttgaaacttgtctatgaggtaggtcctgaaagtttcttgccaccgccaactgtcaaatcagccctgttaaacattaaaagaaaacacttattttttgattttaagtttaaagccaaatacttagcatttatttcctgtctgttagagaaacctgatttatctgtaaaaacagctttaaagtcgattttcaggaaaagtcaggtcaggtcaatttcggaaaaattcggtttaaaccttaatgctcaaattgtttgtttgtctccaagtcaatggttaaactgttttttggaaatgctggaagttgtccctgaaaaatttcatccttcgtagatgacaaaaaagaaattgcccgttcgttttacgggtcagcactttactattgataaagtgctaataaaagatgcaataagacaagcaaatataagtaatcaggatacgg")),
  tar_target(
    name =ln_vsa,
    command = desman_input_file|>get_data_from_server()|>read_desman_input_files()),
  tar_target(
    name=tau_vga,
    command =variants|>
      translate_dna_string_vector_to_string_matrix()|>
      translate_dna_matrix_to_binary_array()),
  tar_target(name = mcmc_output,
             command = desman_fixed_variants(ln_vsa[[1]],
                                                       tau_vga,
                                                       G=dim(tau_vga)[2]-1,
                                             tildeepsilon=NA,
                                                       error_rate = 0.001,
                                                       prior_std = 0.01,
                                                       n_chains=2,
                                                       alpha0=.1)),
  tar_target(name = mcmc_output_fixed_epsilon,
             command = desman_fixed_variants(ln_vsa[[1]],
                                             tau_vga,
                                             G=dim(tau_vga)[2]-1,
                                             tildeepsilon=.0001,
                                             error_rate = 0.001,
                                             prior_std = 0.01,
                                             n_chains=2,
                                             alpha0=.1)))
