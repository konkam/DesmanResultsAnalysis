library(targets)
Sys.setenv(TAR_PROJECT = "fixed_variants")
tar_make()

Sys.setenv(TAR_PROJECT = "fixed_variants_large_dictionnary")
tar_make()
