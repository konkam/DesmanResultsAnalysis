compute_relative_abundance <- function(number_of_variants = 1, chain = 5, include_warmup = FALSE, prefix = "", quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) {
  load_gamma_parameter(number_of_variants = number_of_variants, prefix = prefix, include_warmup = include_warmup, nchains = chain) %>%
    ggmcmc::ggs() %>%
    # dplyr::filter(Chain == chain) %>%
    (function(df) {
      df[df$Chain == chain, ]
    }) %>%
    calculate_quantiles_by_factor(factor_column = "Parameter", value_column = "value", quantiles = quantiles, add_mean = T) %>%
    dplyr::mutate(
      Variant = Parameter %>% stringr::str_split_i("_", i = 2),
      Sample_id = Parameter %>% stringr::str_split_i("_", i = 3)
    )
}

# compute_relative_abundance(number_of_variants = 2, nchains = 5, include_warmup = FALSE, prefix = "test_1/")

#' Plot the relative abundance for each variants in each sample
#'
#' @param number_of_variants The number of variants considered, this is used to choose the file path
#' @param chain The chain selected. This ought to be the chain with the lowest deviance among all chains.
#' @param prefix File path describing where the results are stored.
#' @importFrom dplyr mutate
#' @importFrom ggplot2 geom_segment scale_x_continuous geom_point
#'
#' @return A plot the relative abundance for each variants in each sample
#' @export
#'
plot_relative_abundance <- function(number_of_variants = 1, chain = 5, prefix = "") {
  compute_relative_abundance(number_of_variants = number_of_variants, chain = chain, include_warmup = FALSE, prefix = prefix, quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
    mutate(
      Variant = as.numeric(Variant),
      Sample_id = as.numeric(Sample_id)
    ) %>%
    ggplot(aes(x = Variant)) +
    theme_bw() +
    facet_wrap(~Sample_id) +
    geom_segment(aes(xend = Variant, y = percentile_2.5, yend = percentile_97.5)) +
    geom_segment(aes(xend = Variant, y = percentile_25, yend = percentile_75), alpha = 0.9, colour = "grey", linewidth = 2) +
    geom_point(aes(y = mean)) +
    ylim(0, 1) +
    ylab("Relative abundance") +
    scale_x_continuous(breaks = function(x) seq(ceiling(x[1]), floor(x[2]), by = 1))
}
