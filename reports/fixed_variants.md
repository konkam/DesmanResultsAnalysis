MCMC report
================

    ## [1] "runjags"

    ## ✖ Identified no targets to load.

# Setup

## Project ID

smc_1

## Input data

The counts per position, sample and nucleotides array was obtained from:
[/work_projet/ala/metachick-fugace/Analyses_Ouléye/DESMAN_2023/DESMAN_iterAddPost/DESMAN_erm_F\_\_3_M17808/freq2Desman.txt]()
An ad hoc dictionnary of size 4 was used

## Desman Tuning parameters

10000
1000
1000
0.001
0.01
2
0.1

# Results

## Free epsilon

``` r
plot_tilde_epsilon(mcmc_output = mcmc_output,error_rate = desman_tuning_parameters$error_rate)
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](fixed_variants_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
plot_pi(mcmc_output,n_vsa =n_vsa,variants = variants)
```

![](fixed_variants_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->

## Fixed epsilon

``` r
plot_pi(mcmc_output=mcmc_output_fixed_epsilon,n_vsa =n_vsa,variants = variants )
```

![](fixed_variants_files/figure-gfm/unnamed-chunk-5-1.png)<!-- --> \#
Reproduction To reproduce the code: If not done, clone the
[https://github.com/konkam/DesmanResultsAnalysis](repo)

then run:

``` r
library(targets)
Sys.setenv(TAR_PROJECT = 'smc_1')
tar_make()
```
