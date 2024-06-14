Problem with tau sampler
================

# Simulate data

``` r
tau_vgb_s=c("aaaaaa","aaaccc","cccaaa","cccttt","cccggg")|>
  translate_dna_string_vector_to_string_matrix()
tau_vgb=tau_vgb_s|>
  translate_dna_matrix_to_binary_array()

tau_vgb1_s=c("accaaa","caaggg","caaatt","caataa","accccc")|>
  translate_dna_string_vector_to_string_matrix()
tau_vgb1=tau_vgb1_s|>
  translate_dna_matrix_to_binary_array()

pi_gs<-matrix(c(.2,.3,.1,.1,.2),5,1)

tau_vgb_s
```

    ##    g
    ## v   1   2   3   4   5  
    ##   1 "a" "a" "c" "c" "c"
    ##   2 "a" "a" "c" "c" "c"
    ##   3 "a" "a" "c" "c" "c"
    ##   4 "a" "c" "a" "t" "g"
    ##   5 "a" "c" "a" "t" "g"
    ##   6 "a" "c" "a" "t" "g"

``` r
tau_vgb1_s
```

    ##    g
    ## v   1   2   3   4   5  
    ##   1 "a" "c" "c" "c" "a"
    ##   2 "c" "a" "a" "a" "c"
    ##   3 "c" "a" "a" "a" "c"
    ##   4 "a" "g" "a" "t" "c"
    ##   5 "a" "g" "t" "a" "c"
    ##   6 "a" "g" "t" "a" "c"

``` r
pi_gs
```

    ##      [,1]
    ## [1,]  0.2
    ## [2,]  0.3
    ## [3,]  0.1
    ## [4,]  0.1
    ## [5,]  0.2

``` r
error_rate <- .0001
epsilon_ba <- epsilon_ba_f(error_rate)


n_vsa=sim_n_vsa(tau_vgb=tau_vgb,n = 10000,pi_gs = pi_gs,error_rate = error_rate)
drop(n_vsa)|> (function(x){dimnames(x)[[2]]<-nucleotides;x})()
```

    ##    a
    ## v      a    c    g    t
    ##   1 5558 4441    0    1
    ##   2 5598 4402    0    0
    ##   3 5533 4465    1    1
    ##   4 3369 3269 2234 1128
    ##   5 3251 3390 2223 1136
    ##   6 3268 3333 2253 1146

``` r
g=dim(tau_vgb)[2]
s=1
```

What if we start from the dictionnary that was used to simulate the data
? We do not move:

``` r
plyr::raply(10,sampler_tau(tau_vgb,pi_gs,epsilon_ba,n_vsa,
                      v=dim(tau_vgb)[1],
                      g=dim(pi_gs)[1],
                      g_neq_g=g_neq_g_f(g))|>unname()|>
              translate_dna_binary_array_to_string_vector())|>rbind(
tau_vgb|>
  translate_dna_binary_array_to_string_vector()|>unname())
```

    ##       1        2        3        4        5       
    ##  [1,] "aaaaaa" "aaaccc" "cccaaa" "cccttt" "cccggg"
    ##  [2,] "aaaaaa" "aaaccc" "cccaaa" "cccttt" "cccggg"
    ##  [3,] "aaaaaa" "aaaccc" "cccaaa" "cccttt" "cccggg"
    ##  [4,] "aaaaaa" "aaaccc" "cccaaa" "cccttt" "cccggg"
    ##  [5,] "aaaaaa" "aaaccc" "cccaaa" "cccttt" "cccggg"
    ##  [6,] "aaaaaa" "aaaccc" "cccaaa" "cccttt" "cccggg"
    ##  [7,] "aaaaaa" "aaaccc" "cccaaa" "cccttt" "cccggg"
    ##  [8,] "aaaaaa" "aaaccc" "cccaaa" "cccttt" "cccggg"
    ##  [9,] "aaaaaa" "aaaccc" "cccaaa" "cccttt" "cccggg"
    ## [10,] "aaaaaa" "aaaccc" "cccaaa" "cccttt" "cccggg"
    ## [11,] "aaaaaa" "aaaccc" "cccaaa" "cccttt" "cccggg"

What if we start from a dictionnary that was different from the one used
to simulate the data ? We do not move either:

``` r
plyr::raply(10,sampler_tau(tau_vgb1,pi_gs,epsilon_ba,n_vsa,
                           v=dim(tau_vgb)[1],
                           g=dim(pi_gs)[1],
                           g_neq_g=g_neq_g_f(g))|>unname()|>
              translate_dna_binary_array_to_string_vector())|>rbind(
                tau_vgb1|>
                  translate_dna_binary_array_to_string_vector()|>unname())
```

    ##       1        2        3        4        5       
    ##  [1,] "accaaa" "caaggg" "aaaatt" "aaatcc" "accccc"
    ##  [2,] "accaaa" "caaggg" "aaaatt" "aaatcc" "accccc"
    ##  [3,] "accaaa" "caaggg" "aaaatt" "aaatcc" "accccc"
    ##  [4,] "accaaa" "caaggg" "aaaatt" "aaatcc" "accccc"
    ##  [5,] "accaaa" "caaggg" "aaaatt" "aaatcc" "accccc"
    ##  [6,] "accaaa" "caaggg" "aaaatt" "aaatcc" "accccc"
    ##  [7,] "accaaa" "caaggg" "aaaatt" "aaatcc" "accccc"
    ##  [8,] "accaaa" "caaggg" "aaaatt" "aaatcc" "accccc"
    ##  [9,] "accaaa" "caaggg" "aaaatt" "aaatcc" "accccc"
    ## [10,] "accaaa" "caaggg" "aaaatt" "aaatcc" "accccc"
    ## [11,] "accaaa" "caaggg" "caaatt" "caataa" "accccc"

Simpler simulation

``` r
#################################################################################

tau_vgb=c("ac","ca")|>
  translate_dna_string_vector_to_string_matrix()|>
  translate_dna_matrix_to_binary_array()


rdirichlet(c(100,100),1)
```

    ##           [,1]      [,2]
    ## [1,] 0.4869599 0.5130401

``` r
#################################################################################

tau_vgb=c("aa","cc","ac","ca")|>
  translate_dna_string_vector_to_string_matrix()|>
  translate_dna_matrix_to_binary_array()

pi_gs<-matrix(c(.4,.4,.1,.1),4,1)
pi_gs0<-matrix(c(.1,.1,.4,.4),4,1)



alpha_g=.0001
n_vsa=sim_n_vsa(tau_vgb=tau_vgb,n = 10000,pi_gs = pi_gs,error_rate = .0001)

  xi_vsabg=sampler_xi(tau_vgb=tau_vgb,
                      pi_gs=pi_gs,
                      epsilon_ba = epsilon_ba)
  nu_vsab<-nu_from_xi(xi=xi_vsabg)
  mu_vsab<-mu_from_xi(xi=xi_vsabg)
  pi_gs1<-sampler_pi(mu_vsag = mu_vsab,alpha_g = alpha_g)
```

# Simulate data

``` r
set.seed(1)
tau_pi_n <- sim_tau_pi_epsilon_n(v = 6, g = 5, s = 3, n = 1000, alpha_pi = 1)

tau_pi_n$tau_vgb|> translate_dna_binary_array_to_string_vector()
```

    ##        1        2        3        4        5 
    ## "atgaca" "ggccgg" "aaaccc" "cgagaa" "aacaac"

``` r
tau_pi_n$pi_gs
```

    ##            [,1]       [,2]       [,3]
    ## [1,] 0.18085786 0.02272180 0.23993954
    ## [2,] 0.38898245 0.45343579 0.23162256
    ## [3,] 0.21380266 0.40214570 0.30009026
    ## [4,] 0.05045933 0.07298370 0.06668715
    ## [5,] 0.16589770 0.04871302 0.16166049

``` r
tau_pi_n$n_vsa|>display_n_vsa()
```

    ##    s-a
    ## v   s1-a s1-c s1-g s1-t s2-a s2-c s2-g s2-t s3-a s3-c s3-g s3-t
    ##   1  560   51  387    2  469   76  455    0  699   62  239    0
    ##   2  377    0  446  177  458    0  525   17  445    0  296  259
    ##   3  253  564  183    0  475  490   35    0  381  379  240    0
    ##   4  355  605   39    1   82  841   77    0  423  508   69    0
    ##   5  223  384  391    2  123  433  444    0  248  512  240    0
    ##   6  251  362  386    1   86  465  447    2  290  470  239    1

# Test samplers

## Desman

``` r
 gs="jags"
 block_tau=FALSE
 shape_rho=c(1,100)
 GS=mcmc_very_relax_run(
   n_vsa = tau_pi_n$n_vsa,
   G = 12
 )
```

    ## Warning: No initial value blocks found and n.chains not specified: 2 chains
    ## were used

    ## Le chargement a nécessité le package : rjags

    ## Warning: No initial values were provided - JAGS will use the same initial
    ## values for all chains

    ## Compiling rjags model...
    ## Calling the simulation using the rjags method...
    ## Adapting the model for 1000 iterations...

    ## Warning: The adaptation phase of the model was not completed in 1000
    ## iterations, so the current samples may not be optimal - try increasing the
    ## number of iterations to the "adapt" argument

    ## Burning in the model for 4000 iterations...
    ## Running the model for 10000 iterations...

    ## Warning in FUN(X[[i]], ...): Failed to set trace monitor for tau
    ## Variable tau not found

    ## Simulation complete
    ## Note: Summary statistics were not produced as there are >50 monitored
    ## variables
    ## [To override this behaviour see ?add.summary and ?runjags.options]
    ## FALSEFinished running the simulation

``` r
 GS$mcmc[[1]]->X;
 X[10000,paste0("pi_gs[",c(outer(1:12,1:3,paste,sep=",")),"]")]
```

    ##  pi_gs[1,1]  pi_gs[2,1]  pi_gs[3,1]  pi_gs[4,1]  pi_gs[5,1]  pi_gs[6,1] 
    ## 0.008123673 0.118623299 0.071000344 0.068677835 0.180382323 0.016559739 
    ##  pi_gs[7,1]  pi_gs[8,1]  pi_gs[9,1] pi_gs[10,1] pi_gs[11,1] pi_gs[12,1] 
    ## 0.134835957 0.028959262 0.053874812 0.092165749 0.164541462 0.062255545 
    ##  pi_gs[1,2]  pi_gs[2,2]  pi_gs[3,2]  pi_gs[4,2]  pi_gs[5,2]  pi_gs[6,2] 
    ## 0.048768366 0.020479927 0.004518537 0.028862091 0.020092973 0.000917146 
    ##  pi_gs[7,2]  pi_gs[8,2]  pi_gs[9,2] pi_gs[10,2] pi_gs[11,2] pi_gs[12,2] 
    ## 0.030317953 0.068488541 0.038380664 0.047353984 0.326637489 0.365182329 
    ##  pi_gs[1,3]  pi_gs[2,3]  pi_gs[3,3]  pi_gs[4,3]  pi_gs[5,3]  pi_gs[6,3] 
    ## 0.048020910 0.153441654 0.073633461 0.046713577 0.218687772 0.022227358 
    ##  pi_gs[7,3]  pi_gs[8,3]  pi_gs[9,3] pi_gs[10,3] pi_gs[11,3] pi_gs[12,3] 
    ## 0.172246062 0.019344509 0.013931934 0.028046559 0.013842667 0.189863538

``` r
 A=matrix(NA,12,3) ;for (i in 1:12){for(j in 1:3){A[i,j]=X[10000,paste0("pi_gs[",i,",",j,"]")]}}
```
