collect_variants_one_chain <- function(tau_sample_one_chain) {
  # tau_sample_one_chain = load_tau_one_chain(number_of_variants = 3, chain = 1, include_warmup = T, prefix = path_desman_results)
  # The nucleotide maximising the L2 posterior loss is the most frequent nucleotide, or the MAP of the MCMC sample.
  n_variants <- dim(tau_sample_one_chain)[2]
  n_positions <- dim(tau_sample_one_chain)[3]

  tau_sample_one_chain %>%
    apply(X = ., MARGIN = c(1, 2, 3), FUN = sum) %>%
    apply(X = ., MARGIN = c(2, 3), FUN = which.max) %>%
    apply(X = ., c(1, 2), function(x) nucleotides_letters[x]) %>%
    t() %>%
    as_tibble() %>%
    rowid_to_column("Position")
}

merge_variants_all_chains <- function(tau_samples_all_chains, keep_only_variable_positions = F) {
  # tau_samples_all_chains = lapply(1:10, function(chain_id){load_tau_one_chain(number_of_variants = 3, chain = chain_id, include_warmup = T, prefix = path_desman_results)})

  # lapply(seq_along(tau_samples_all_chains), function(chain_id){
  #   collect_variants_one_chain(tau_samples_all_chains[[chain_id]]) %>%
  #     mutate(chain_id = chain_id)
  # })

  # tau_samples_all_chains %>%
  #   lapply(collect_variants_one_chain) %>%
  #   lapply(FUN = function(df) df %>% select(-Position)) %>%
  #   bind_cols() %>%
  #   select(across(everything(), ~!duplicated(.)))
  #
  # tau_samples_all_chains %>%
  #   lapply(collect_variants_one_chain) %>%
  #   lapply(FUN = function(df) df %>% select(-Position)) %>%
  #   bind_cols() %>%
  #   mutate(across(everything(), ~ ifelse(duplicated(.), NA, .))) %>%
  #   select(everything(), -where(~all(is.na(.))))

  res <- tau_samples_all_chains %>%
    lapply(collect_variants_one_chain) %>%
    lapply(FUN = function(df) df %>% select(-Position)) %>%
    bind_cols() %>%
    t() %>%
    as_tibble() %>%
    distinct() %>%
    t() %>%
    as_tibble() %>%
    rowid_to_column("Position")

  if (keep_only_variable_positions) {
    res %>%
      filter_rows_not_all_equal(columns_to_exclude = "Position")
  } else {
    res
  }
}

which_variants_in_each_chain <- function(tau_samples_all_chains) {
  reference_all_positions <- merge_variants_all_chains(tau_samples_all_chains, keep_only_variable_positions = F)
  reference_variable_positions <- merge_variants_all_chains(tau_samples_all_chains, keep_only_variable_positions = T)
  seq_along(tau_samples_all_chains) %>%
    lapply(
      (function(chain_id) {
        variants_in_chain <- collect_variants_one_chain(tau_samples_all_chains[[chain_id]])
        n_vars <- variants_in_chain %>%
          ncol() %>%
          (function(x) x - 1)
        lapply(
          1:n_vars,
          (function(variant_id_in_chain) {
            tibble(chain_id = chain_id, variant_id_in_chain = variant_id_in_chain, variant_id_in_ref = find_equal_column(reference_all_positions, variants_in_chain[, paste("V", variant_id_in_chain, sep = "")]) %>%
              gsub("V", "", .) %>%
              as.numeric())
          })
        ) %>%
          bind_rows()
      })
    ) %>%
    bind_rows() %>%
    mutate(
      Variant_full_sequence = lapply(variant_id_in_ref, function(vid) {
        reference_all_positions[, paste("V", vid, sep = "")] %>%
          unlist() %>%
          Reduce(function(ncl1, ncl2) paste(ncl1, ncl2, sep = ""), .)
      }) %>% unlist(),
      Variant_sequence_variable_positions = lapply(variant_id_in_ref, function(vid) {
        reference_variable_positions[, paste("V", vid, sep = "")] %>%
          unlist() %>%
          Reduce(function(ncl1, ncl2) paste(ncl1, ncl2, sep = ""), .)
      }) %>% unlist()
    )
}
