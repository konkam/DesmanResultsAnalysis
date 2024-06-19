
#' @examples
#' dna_string <- "ActGG"
#' result <- translate_dna_to_binary(dna_string)
#' print(result)
translate_dna_to_binary <- function(dna_string) {
  # Create a matrix to store the binary representation
  binary_matrix <- matrix(0,
                          nrow = dna_string|>nchar(),
                          ncol = 4,
                          dimnames = list(NULL, nucleotides)
  )
  
  # Fill in the matrix based on the nucleotides in the DNA string
  for (i in 1:nchar(dna_string)) {
    nucleotide <- dna_string|>tolower()|>substr(i, i)
    binary_matrix[i, nucleotide] <- 1
  }
  
  return(binary_matrix)
}


#' @examples
#' variants_string_matrix <- sim_tau_vgb_string_matrix(3, 12)
#' translate_dna_matrix_to_binary_array(variants_string_matrix)
translate_dna_matrix_to_binary_array <- function(variants_string_matrix) {
  variants_string_matrix |>
    tolower() |>
    plyr::aaply(1:2, `==`, nucleotides) |>
    (`*`)(1)|>
    (function(x){names(dimnames(x))[3]="b";x})()
}

#' @examples
#' variants_string_vector <- 
#' sim_tau_vgb_string_matrix(g = 3, v = 12) |> 
#' plyr::aaply(2, paste0, collapse = "")
#' translate_dna_string_vector_to_string_matrix(variants_string_vector)
translate_dna_string_vector_to_string_matrix <- function(variants_string_vector) {
  variants_string_vector |>
    tolower() |>
    plyr::aaply(1, function(x) {
      strsplit(x, "") |>
        unname() |>
        unlist()
    }) |>
    t() |>
    (function(x) {
      x |> (`dimnames<-`)(list(v = 1:dim(x)[1], 
                               g = 1:dim(x)[2]))
    })()
}


#' @examples
#' variants_string_vector <- 
#' sim_tau_vgb_string_matrix(g = 3, v = 12) |> 
#' plyr::aaply(2, paste0, collapse = "")
#' translate_dna_string_vector_to_string_matrix(variants_string_vector)|>
#' translate_dna_matrix_to_binary_array()|>
#' translate_dna_binary_array_to_string_vector()
translate_dna_binary_array_to_string_vector <- 
  function(variants_binary_array) {
    variants_binary_array |>
      plyr::aaply(1:2,function(x){nucleotides[x==1]},.drop=FALSE)|>
      abind::adrop(3)|>
      plyr::aaply(2,paste,collapse="")
  }