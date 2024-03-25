#' @example
#' "/work_projet/ala/metachick-fugace/BD/allTEM.fasta"|>
#' get_data_from_server()|>
#' read_fasta_file()->dic
#' dic|>nchar()|>table()
read_fasta_file <- function(path_to_fasta_file) {
  path_to_fasta_file |>
    readLines() |>
    (function(x) {
      start <- grep(x = x, pattern = ">")
      the_names <- x[start] |>
        gsub(pattern = ">", replacement = "") |>
        sub(pattern = "_[^_]+$", replacement = "")
      1:length(start) |>
        plyr::aaply(1, function(i) {
          x[(start[i] + 1):(c(start[-1], (length(x) + 1))[i] - 1)] |>
            paste(collapse = "")
        }) |>
        setNames(the_names)
    })()
}
