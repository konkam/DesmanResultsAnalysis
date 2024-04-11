MCMC report
================

    ## [1] "runjags"

# Setup

## Project ID

fixed_variants_large_dictionnary

## Input data

The counts per position, sample and nucleotides array was obtained from:
[/work_projet/ala/metachick-fugace/Analyses_Ouléye/DESMAN_2023/DESMAN_iterAddPost/DESMAN_blaTEM_1A_1_HM749966/freq2Desman.txt]()
An ad hoc dictionnary of size 176 was used

## MCMC Tuning parameters

# Results

``` r
#plot_tilde_epsilon(mcmc_output = mcmc_output)
```

# Reproduction

To reproduce the code: If not done, clone the
[https://github.com/konkam/DesmanResultsAnalysis](repo)

then run:

``` r
library(targets)
Sys.setenv(TAR_PROJECT = 'fixed_variants_large_dictionnary')
tar_make()
```
