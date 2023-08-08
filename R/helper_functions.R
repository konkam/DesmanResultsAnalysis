#' Calculate quantiles by factor
#'
#' This function calculates quantiles for a given numeric column within each factor group in the input table. The result is a data frame with quantile values in separate columns for each factor.
#'
#' @param input_table A data frame containing the input data with factors and numeric values.
#' @param factor_column A character string specifying the name of the column containing factors.
#' @param value_column A character string specifying the name of the column containing numeric values.
#' @param quantiles A numeric vector of quantiles to compute, e.g., c(0.25, 0.5, 0.75).
#' @return A data frame with quantile values in separate columns for each factor.
#' @examples
#' input_table <- data.frame(
#'   factor_column = c("A", "A", "B", "B", "A", "B"),
#'   value_column = c(10, 20, 15, 25, 30, 35)
#' )
#' quantiles_to_compute <- c(0.25, 0.5, 0.75)
#' result <- calculate_quantiles_by_factor(input_table, "factor_column", "value_column", quantiles_to_compute)
#' print(result)
#'
#' @export
calculate_quantiles_by_factor <- function(input_table, factor_column, value_column, quantiles, add_mean = F) {
  # Ensure that 'quantiles' is sorted in ascending order
  quantiles <- sort(quantiles)

  # Split the data by the factor_column
  split_data <- split(input_table[[value_column]], input_table[[factor_column]])

  # Calculate quantiles for each factor
  quantile_data <- lapply(split_data, function(x) quantile(x, probs = quantiles))

  # Combine the quantile results into a data frame
  output_table <- as.data.frame(t(sapply(quantile_data, unlist)))
  colnames(output_table) <- paste0("percentile_", quantiles * 100)

  output_table = output_table %>%
    (function(df){
      df %>%
        as_tibble() %>%
        bind_cols(tibble(!!sym(factor_column)  := df %>% rownames()), .)
    })

  if(add_mean){
    mean_table = split_data %>%
      lapply(mean) %>%
      unlist %>%
      (function(l){
        tibble(!!sym(factor_column)  := names(l), mean = l)
      })

    output_table = left_join(output_table, mean_table)
  }

  return(output_table)
}
