#' Calculate quantiles by factor
#'
#' This function calculates quantiles for a given numeric column within each factor group in the input table. The result is a data frame with quantile values in separate columns for each factor.
#'
#' @param input_table A data frame containing the input data with factors and numeric values.
#' @param factor_column A character string specifying the name of the column containing factors.
#' @param value_column A character string specifying the name of the column containing numeric values.
#' @param quantiles A numeric vector of quantiles to compute, e.g., c(0.25, 0.5, 0.75).
#' @return A data frame with quantile values in separate columns for each factor.
#' @importFrom tibble tibble
#' @importFrom dplyr sym
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
calculate_quantiles_by_factor <- function(input_table, factor_column, value_column, quantiles, add_mean = FALSE) {
  # Ensure that 'quantiles' is sorted in ascending order
  quantiles <- sort(quantiles)

  # Split the data by the factor_column
  split_data <- split(input_table[[value_column]], input_table[[factor_column]])

  # Calculate quantiles for each factor
  quantile_data <- lapply(split_data, function(x) quantile(x, probs = quantiles))

  # Combine the quantile results into a data frame
  output_table <- as.data.frame(t(sapply(quantile_data, unlist)))
  colnames(output_table) <- paste0("percentile_", quantiles * 100)

  output_table <- output_table %>%
    (function(df) {
      df %>%
        tibble::as_tibble() %>%
        dplyr::bind_cols(tibble(!!sym(factor_column) := df %>% rownames()), .)
    })

  if (add_mean) {
    mean_table <- split_data %>%
      lapply(mean) %>%
      unlist() %>%
      (function(l) {
        tibble(!!sym(factor_column) := names(l), mean = l)
      })

    output_table <- dplyr::left_join(output_table, mean_table)
  }

  return(output_table)
}

nucleotides_letters <- c("A", "C", "G", "T")
nucleotides_indicator_vector = lapply(1:4, function(i){diag(x = 1, nrow = 4)[i,]})
nucleotides_dict = c(nucleotides_letters, nucleotides_indicator_vector) %>% (function(v) v %>% setNames(rev(v)))
nucleotides_conversion_function = function(nms){
  sapply(nms, function(nm){nucleotides_dict[nm]})
}


#' Title
#'
#' @param input_tibble
#' @param columns_to_check
#'
#' @return
#'
#' @examples
#' # Sample data
#' df <- tibble(
#'   ID = c(1, 2, 3, 4, 5, 6),
#'   A = c(10, 20, 10, 30, 5, 6),
#'   B = c(5, 5, 5, 5, 5, 5),
#'   C = c("X", "Y", "X", "Z", "Z", "Z")
#' )
#'
#' #Only check columns A and B
#' result_df1 <- filter_rows_not_all_equal(df, columns_to_check = c("A", "B"))
#' print(result_df1)
#'
#' # Check all columns except ID
#' result_df2 <- filter_rows_not_all_equal(df, columns_to_exclude = c("C", "ID"))
#' print(result_df2)
#'
#' # Check all columns
#' result_df3 <- filter_rows_not_all_equal(df)
#' print(result_df3)
filter_rows_not_all_equal <- function(input_tibble, columns_to_check = NULL, columns_to_exclude = NULL) {
  if (!is.null(columns_to_check)) {
    result_tibble <- input_tibble %>%
      filter(
        (!do.call(pmin, dplyr::select(., dplyr::all_of(columns_to_check))) ==
           do.call(pmax, dplyr::select(., dplyr::all_of(columns_to_check))))
        | is.na(do.call(pmin, dplyr::select(., dplyr::all_of(columns_to_check))))
      )
  } else {
    result_tibble <- input_tibble %>%
      filter(
        (!do.call(pmin, dplyr::select(., -dplyr::all_of(columns_to_exclude))) ==
           do.call(pmax, dplyr::select(., -dplyr::all_of(columns_to_exclude))))
        | is.na(do.call(pmin, dplyr::select(., -dplyr::all_of(columns_to_exclude))))
      )
  }
  return(result_tibble)
}

#' Find column name equal to vector, if ant
#'
#' @param input_tibble
#' @param target_vector
#'
#' @return A column name, or character(0) if there is no match
#' @export
#'
#' @examples
#' # Sample data
#' df <- tibble(
#'   A = c(1, 2, 1),
#'   B = c(4, 5, 6),
#'   C = c(7, 8, 9)
#' )
#'
#' # Target vector to compare
#' target_vector <- c(1, 2, 1)
#'
#' # Find columns equal to the target vector
#' matching_columns <- find_equal_column(df, target_vector)
#'
#' # Print the result
#' print(matching_columns)
find_equal_column <- function(input_tibble, target_vector) {
  matching_columns <- colnames(input_tibble)[sapply(input_tibble, function(col) all(col == target_vector))]
  return(matching_columns)
}

