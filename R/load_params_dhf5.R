# library(coda)
# library(tidyverse)

create_folder_name <- function(number_of_variants, chain, prefix = "") {
  paste(prefix, "desman.", number_of_variants, ".", chain, sep = "")
}

load_eta_parameter <- function(number_of_variants = 1, nchains = 5, include_warmup = FALSE, prefix = "") {
  1:nchains %>%
    lapply(FUN = function(chain_id) {
      load_eta_one_chain(number_of_variants = number_of_variants, chain = chain_id, include_warmup = include_warmup, prefix = prefix)
    }) %>%
    coda::mcmc.list()
}

load_eta_one_chain <- function(number_of_variants = 1, chain = 1, include_warmup = FALSE, prefix = "") {
  foldername <- create_folder_name(number_of_variants = number_of_variants, chain = chain, prefix = prefix)

  # warm_up_eta <- hdf5r::H5File$new(paste(foldername, "/eta_store_before_burnin.csv", sep = ""), mode = "r+") %>%
  #   (function(fileh5) {
  #     fileh5[["eta_store"]][, , ]
  #   })
  #
  # sampling_eta <- hdf5r::H5File$new(paste(foldername, "/eta_store.csv", sep = ""), mode = "r+") %>%
  #   (function(fileh5) {
  #     fileh5[["eta_store"]][, , ]
  #   })

  warm_up_eta <- rhdf5::h5read(paste(foldername, "/eta_store_before_burnin.csv", sep = ""), name = "eta_store")
  sampling_eta <- rhdf5::h5read(paste(foldername, "/eta_store.csv", sep = ""), name = "eta_store")

  mcmc_object <- tidyr::expand_grid(i = 1:4, j = 1:4) %>%
    (function(df) {
      dplyr::bind_rows(
        mapply(
          FUN = function(i, j) {
            tibble(eta = warm_up_eta[i, j, ]) %>%
              setNames(paste("eta_", i, "_", j, sep = ""))
          },
          df$i, df$j, SIMPLIFY = FALSE
        ) %>%
          dplyr::bind_cols(), # %>%
        # rowid_to_column(var = "iter")
        mapply(
          FUN = function(i, j) {
            tibble(eta = sampling_eta[i, j, ]) %>%
              setNames(paste("eta_", i, "_", j, sep = ""))
          },
          df$i, df$j, SIMPLIFY = FALSE
        ) %>%
          dplyr::bind_cols()
      )
    }) %>%
    # select(-iter) %>%
    (function(df) {
      if (!include_warmup) {
        df[round(nrow(df) / 2 + 1):nrow(df), ]
      } else {
        df
      }
    }) %>%
    coda::mcmc(data = .)

  return(mcmc_object)
}

load_ll <- function(number_of_variants = 1, nchains = 5, include_warmup = FALSE, prefix = "") {
  1:nchains %>%
    lapply(FUN = function(chain_id) {
      foldername <- create_folder_name(number_of_variants, chain_id, prefix = prefix)
      res <- readr::read_csv(
        file = paste(foldername, "/ll_store.csv", sep = ""),
        col_names = c("iter", "loglik")
      ) %>%
        .[-1, "loglik"]
      if (include_warmup) {
        dplyr::bind_rows(
          readr::read_csv(
            file = paste(foldername, "/ll_store_before_burnin.csv", sep = ""),
            col_names = c("iter", "loglik")
          ) %>%
            .[-1, "loglik"],
          res
        ) %>%
          coda::mcmc()
      } else {
        res %>% coda::mcmc()
      }
    }) %>%
    coda::mcmc.list()
}

nucleotides_letters <- c("A", "C", "G", "T")

summarise_tau <- function(tau) {
  # if (length(dim(tau)) == 3) {
  #   number_of_variants <- 1
  # } else {
  #   number_of_variants <- dim(tau)[2]
  # }

  number_of_variants <- dim(tau)[2]

  tidyr::expand_grid(i = seq_along(nucleotides_letters), j = 1:number_of_variants) %>%
    (function(df) {
      mapply(FUN = function(nucleotide_id, variant_id) {
        tmp <- tau[nucleotide_id, variant_id, , ] %>% colSums() # Counting the number of nucleotide of type nucleotide_id in the genome of variant_id, at each iteration.
        tibble(aa = tmp) %>%
          setNames(paste(nucleotides_letters[nucleotide_id], variant_id, sep = ""))
      }, df$i, df$j, SIMPLIFY = FALSE)
    }) %>%
    dplyr::bind_cols()
}

load_tau_one_chain <- function(number_of_variants = 1, chain = 1, include_warmup = FALSE, prefix = "") {
  foldername <- create_folder_name(number_of_variants, chain, prefix = prefix)

  # Note that if you have 1 variant a 3d array is loaded, while with more than 1 variants it's a 4d variable
  # 4 X number_of_variants X number of variable positions X number of iterations
  # if (include_warmup) {
  #   warm_up_tau <- hdf5r::H5File$new(paste(foldername, "/tau_store_before_burnin.csv", sep = ""), mode = "r+") %>%
  #     (function(fileh5) {
  #       fileh5[["tau_store"]][, , , ]
  #     })
  # }
  #
  # sampling_tau <- hdf5r::H5File$new(paste(foldername, "/tau_store.csv", sep = ""), mode = "r+") %>%
  #   (function(fileh5) {
  #     fileh5[["tau_store"]][, , , ]
  #   })

  if (include_warmup) {
    warm_up_tau <- rhdf5::h5read(paste(foldername, "/tau_store_before_burnin.csv", sep = ""), name = "tau_store")
  }

  sampling_tau <- rhdf5::h5read(paste(foldername, "/tau_store.csv", sep = ""), name = "tau_store")

  mcmc_object <- summarise_tau(sampling_tau) %>%
    # select(-iter) %>%
    (function(df) {
      if (include_warmup) {
        dplyr::bind_rows(summarise_tau(warm_up_tau), df) # put warmup in front of course
      } else {
        df
      }
    }) %>%
    coda::mcmc(data = .)

  return(mcmc_object)
}

load_tau_parameter <- function(number_of_variants = 1, nchains = 5, include_warmup = FALSE, prefix = "") {
  1:nchains %>%
    lapply(FUN = function(chain_id) {
      load_tau_one_chain(number_of_variants = number_of_variants, chain = chain_id, include_warmup = include_warmup, prefix = prefix)
    }) %>%
    coda::mcmc.list()
}


load_gamma_one_chain <- function(number_of_variants = 1, chain = 1, include_warmup = FALSE, prefix = "") {
  foldername <- create_folder_name(number_of_variants = number_of_variants, chain = chain, prefix = prefix)

  if (include_warmup) {
    # warm_up_gamma <- hdf5r::H5File$new(paste(foldername, "/Gamma_store_before_burnin.csv", sep = ""), mode = "r+") %>%
    #   (function(fileh5) {
    #     fileh5[["gamma_store"]][, , ]
    #   })

    warm_up_gamma <- rhdf5::h5read(paste(foldername, "/Gamma_store_before_burnin.csv", sep = ""), name = "gamma_store")


  #   if (number_of_variants == 1) {
  #     warm_up_gamma <- warm_up_gamma %>% (function(x) array(data = x, dim = c(1, dim(x)))) # resizing to a tensor for the following code to still work. Note that with 1 variant, proportions are always 1, so not very interesting.
  #   }
  }

  # sampling_gamma <- hdf5r::H5File$new(paste(foldername, "/Gamma_store.csv", sep = ""), mode = "r+") %>%
  #   (function(fileh5) {
  #     fileh5[["gamma_store"]][, , ]
  #   })

  sampling_gamma <- rhdf5::h5read(paste(foldername, "/Gamma_store.csv", sep = ""), name = "gamma_store")


  # if (number_of_variants == 1) {
  #   sampling_gamma <- sampling_gamma %>% (function(x) array(data = x, dim = c(1, dim(x)))) # resizing to a tensor for the following code to still work. Note that with 1 variant, proportions are always 1, so not very interesting.
  # }

  n_samples <- dim(sampling_gamma)[2]


  mcmc_object <- tidyr::expand_grid(i = 1:number_of_variants, j = 1:n_samples) %>%
    (function(df) {
      mapply(
        FUN = function(i, j) {
          res <- tibble(eta = sampling_gamma[i, j, ]) %>%
            setNames(paste("gamma_", i, "_", j, sep = ""))

          if (include_warmup) {
            res <- dplyr::bind_rows(tibble(eta = warm_up_gamma[i, j, ]) %>%
              setNames(paste("gamma_", i, "_", j, sep = "")), res)
          }
          return(res)
        },
        df$i, df$j, SIMPLIFY = FALSE
      ) %>%
        dplyr::bind_cols()
    }) %>%
    coda::mcmc(data = .)

  return(mcmc_object)
}

load_gamma_parameter <- function(number_of_variants = 1, nchains = 5, include_warmup = FALSE, prefix = "") {
  1:nchains %>%
    lapply(FUN = function(chain_id) {
      load_gamma_one_chain(number_of_variants = number_of_variants, chain = chain_id, include_warmup = include_warmup, prefix = prefix)
    }) %>%
    coda::mcmc.list()
}
