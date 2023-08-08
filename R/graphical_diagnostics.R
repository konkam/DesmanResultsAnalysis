# library(MCMCvis)
# library(ggmcmc)
# library(coda)
# source("R/load_params_dhf5.R")

log_likelihood_trace <- function(number_of_variants = 1, nchains = 5, include_warmup = F, prefix = "") {
  load_ll(number_of_variants = number_of_variants, nchains = nchains, include_warmup = include_warmup, prefix = prefix) %>%
    traceplot()
}

eta_trace <- function(number_of_variants = 1, nchains = 5, include_warmup = F, prefix = "") {
  load_eta_parameter(number_of_variants = number_of_variants, nchains = nchains, include_warmup = include_warmup, prefix = prefix) %>%
    ggmcmc::ggs() %>%
    ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
    theme_bw() +
    facet_wrap(~Parameter) +
    geom_line() +
    ylab("Value") +
    viridis::scale_color_viridis(name = "Chain", discrete = T)
}

nucleotide_count_trace = function(number_of_variants = 1, nchains = 5, include_warmup = F, prefix = "") {
  load_tau_parameter(number_of_variants = number_of_variants, nchains = nchains, include_warmup = include_warmup, prefix = prefix) %>%
    ggmcmc::ggs() %>%
    ggplot(aes(x = Iteration, y = value + 0.05*Chain, colour = factor(Chain))) +
    theme_bw() +
    facet_wrap(~Parameter, scales = "free_y") +
    geom_line() +
    ylab("Value") +
    viridis::scale_color_viridis(name = "Chain", discrete = T)
}

gamma_trace <- function(number_of_variants = 1, nchains = 5, include_warmup = F, prefix = "", variants = seq_len(number_of_variants), samples = c(1)) {
  load_gamma_parameter(number_of_variants = number_of_variants, prefix = prefix, include_warmup = include_warmup, nchains = nchains) %>%
    ggmcmc::ggs() %>%
    mutate(Variant = Parameter %>% str_split_i("_", i = 2),
           Sample_id = Parameter %>% str_split_i("_", i = 3)) %>%
    filter(Variant %in% variants) %>%
    filter(Sample_id %in% samples)  %>%
    ggplot(aes(x = Iteration, y = value, colour = factor(Chain))) +
    theme_bw() +
    facet_grid(Variant~Sample_id, labeller = "label_both") +
    geom_line() +
    ylab("Relative abundance") +
    viridis::scale_color_viridis(name = "Chain", discrete = T) +
    ylim(0,1)

}

# calculate_quantiles_by_factor <- function(input_table, factor_column, value_column, quantiles) {
#   # Ensure that 'quantiles' is sorted in ascending order
#   quantiles <- sort(quantiles)
#
#   # Split the data by the factor_column
#   split_data <- split(input_table[[value_column]], input_table[[factor_column]])
#
#   # Calculate quantiles for each factor
#   quantile_data <- lapply(split_data, function(x) quantile(x, probs = quantiles))
#
#   # Combine the quantile results into a data frame
#   output_table <- as.data.frame(t(sapply(quantile_data, unlist)))
#   colnames(output_table) <- paste0("quantile_", round(quantiles * 100))
#
#   return(output_table)
# }
# #
# # # Example usage
# # input_table <- data.frame(
# #   factor_column = c("A", "A", "B", "B", "A", "B"),
# #   value_column = c(10, 20, 15, 25, 30, 35)
# # )
# #
# # quantiles_to_compute <- c(0.25, 0.5, 0.75)
# #
# # result <- calculate_quantiles_by_factor(input_table, "factor_column", "value_column", quantiles_to_compute)
# # print(result)
#
#
# load_gamma_parameter(number_of_variants = 2, prefix = prefix, include_warmup = T, nchains = 5) %>%
#   ggmcmc::ggs() %>%
#   quantile_table(factor_col = Parameter, value_col = value, probs = c(0.2, 0.3))
#
# load_gamma_parameter(number_of_variants = 2, prefix = prefix, include_warmup = T, nchains = 5) %>%
#   ggmcmc::ggs() %>%
#   calculate_quantiles_by_factor(factor_col = "Parameter", value_col = "value", quantiles = c(0.2, 0.3))
