
#' Traceplot for the log likelihood
#'
#' @param number_of_variants The number of variants considered.
#' @param nchains Number of chains to display (the code will display chains 1:nchains)
#' @param include_warmup Whether to display the warm up iterations as well
#' @param prefix The path to the folder results
#'
#' @return A trace plot for the log likelihood
#' @export
log_likelihood_trace <- function(number_of_variants = 1, nchains = 5, include_warmup = FALSE, prefix = "") {
  log_likelihood_ggplot2(number_of_variants = number_of_variants, nchains = nchains, include_warmup = include_warmup, prefix = prefix)
}

log_likelihood_trace_coda <- function(number_of_variants = 1, nchains = 5, include_warmup = FALSE, prefix = "") {
  load_ll(number_of_variants = number_of_variants, nchains = nchains, include_warmup = include_warmup, prefix = prefix) %>%
    coda::traceplot()
}

#' @importFrom ggplot2 ggplot aes theme_bw geom_line ylab xlab
log_likelihood_ggplot2 <- function(number_of_variants = 1, nchains = 5, include_warmup = FALSE, prefix = "") {
  load_ll(number_of_variants = number_of_variants, nchains = nchains, include_warmup = include_warmup, prefix = prefix) %>%
    (function(mcmc_list){
      seq_along(mcmc_list) %>%
        lapply(function(chain_id){
          mcmc_list[[chain_id]] %>%
            tibble::as_tibble() %>%
            mutate(chain_id = chain_id) %>%
            tibble::rowid_to_column(var = "iteration")

        })
    }) %>%
    dplyr::bind_rows() %>%
    ggplot(aes(x = iteration, y = loglik, colour = factor(chain_id))) +
    theme_bw() +
    geom_line() +
    viridis::scale_color_viridis(discrete = T, name = "Chain") +
    xlab("Iteration") +
    ylab("Log likelihood")
}

#' Traceplot for the sequencing error matrix
#'
#' @inheritParams log_likelihood_trace
#' @importFrom ggplot2 ggplot aes theme_bw geom_line ylab
#' @return A trace plot for the sequencing error matrix
#' @export
#'
#' @importFrom ggplot2 facet_wrap
#'
eta_trace <- function(number_of_variants = 1, nchains = 5, include_warmup = FALSE, prefix = "") {
  load_eta_parameter(number_of_variants = number_of_variants, nchains = nchains, include_warmup = include_warmup, prefix = prefix) %>%
    ggmcmc::ggs() %>%
    ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
    theme_bw() +
    facet_wrap(~Parameter) +
    geom_line() +
    ylab("Value") +
    viridis::scale_color_viridis(name = "Chain", discrete = TRUE)
}

#' Trace plot for the count of each type of nucleotide
#'
#' @inheritParams log_likelihood_trace
#'
#' @return A trace plot for the count of each type of nucleotide
#' @export
#'
#'
nucleotide_count_trace <- function(number_of_variants = 1, nchains = 5, include_warmup = FALSE, prefix = "") {
  load_tau_parameter(number_of_variants = number_of_variants, nchains = nchains, include_warmup = include_warmup, prefix = prefix) %>%
    ggmcmc::ggs() %>%
    ggplot(aes(x = Iteration, y = value + 0.05 * Chain, colour = factor(Chain))) +
    theme_bw() +
    facet_wrap(~Parameter, scales = "free_y") +
    geom_line() +
    ylab("Value") +
    viridis::scale_color_viridis(name = "Chain", discrete = TRUE)
}

#' Trace plot for the relative abundances
#'
#' @inheritParams log_likelihood_trace
#' @param variants A integer vector describing which variants to plot. By default, all variants.
#' @param samples A integer vector describing which samples to plot. By default, only the first sample.
#'
#' @return A trace plot for the relative abundances
#' @export
#' @importFrom ggplot2 facet_grid ylim
#'
gamma_trace <- function(number_of_variants = 1, nchains = 5, include_warmup = FALSE, prefix = "", variants = seq_len(number_of_variants), samples = c(1)) {
  load_gamma_parameter(number_of_variants = number_of_variants, prefix = prefix, include_warmup = include_warmup, nchains = nchains) %>%
    ggmcmc::ggs() %>%
    mutate(
      Variant = Parameter %>% stringr::str_split_i("_", i = 2),
      Sample_id = Parameter %>% stringr::str_split_i("_", i = 3)
    ) %>%
    # dplyr::filter( {{ Variant }} %in% variants) %>%
    # dplyr::filter({{ Sample_id }} %in% samples) %>%
    (function(df){
      df[df[["Variant"]] %in% as.character(variants),]
    }) %>%
    (function(df){
      df[df[["Sample_id"]] %in% samples,]
    }) %>%
    ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
    theme_bw() +
    facet_grid(Variant ~ Sample_id, labeller = "label_both") +
    geom_line() +
    ylab("Relative abundance") +
    viridis::scale_color_viridis(name = "Chain", discrete = TRUE) +
    ylim(0, 1)
}

#' Trace plot for the relative abundances for all chains, with chain-specific variants disambiguated.
#'
#' @inheritParams log_likelihood_trace
#' @param variants A integer vector describing which variants to plot. By default, all variants.
#' @param samples A integer vector describing which samples to plot. By default, only the first sample.
#'
#' @return A trace plot for the relative abundances
#' @export
#' @importFrom ggplot2 facet_grid ylim
#'
gamma_trace_with_variants_identified <- function(number_of_variants = 1, nchains = 5, include_warmup = FALSE, prefix = "", variants = NULL, samples = c(1)) {
  variant_table = lapply(1:nchains, function(chain_id){load_tau_one_chain(number_of_variants = number_of_variants, chain = chain_id, include_warmup = include_warmup, prefix = prefix)}) %>%
    which_variants_in_each_chain %>%
    dplyr::rename(Variant = variant_id_in_chain,
           Sample_id = chain_id) %>%
    mutate(Variant = as.character(Variant),
           Sample_id = as.character(Sample_id))

  if(is.null(variants)){
    variants = variant_table$variant_id_in_ref %>% unique()
  }

  load_gamma_parameter(number_of_variants = number_of_variants, prefix = prefix, include_warmup = include_warmup, nchains = nchains) %>%
    ggmcmc::ggs() %>%
    mutate(
      Variant = Parameter %>% stringr::str_split_i("_", i = 2),
      Sample_id = Parameter %>% stringr::str_split_i("_", i = 3)
    ) %>%
    left_join(variant_table) %>%
    (function(df){
      df[df[["variant_id_in_ref"]] %in% as.character(variants),]
    }) %>%
    (function(df){
      df[df[["Sample_id"]] %in% samples,]
    }) %>%
    ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
    theme_bw() +
    facet_grid(Variant_sequence_variable_positions ~ Sample_id) +
    geom_line() +
    ylab("Relative abundance") +
    viridis::scale_color_viridis(name = "Chain", discrete = TRUE) +
    ylim(0, 1)

}
