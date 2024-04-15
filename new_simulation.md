New_simulation
================
2024-04-15

# 1. sparsity 강도가 약한 데이터

precision matrix에 대한 시뮬레이션 데이터를 생성하고, 이에 대한 역행렬을
true precision matrix로 설정하였다.

``` r
N <- 10000
p <- 100
real_pi <- 0.05
sim_data <- simulate_sparse_matrix(real_pi, p)
sim_omega <- round(sim_data, digits = 3)
sim_sigma <- solve(sim_omega)
Y <- mvrnorm(n = N, mu = rep(0, p), Sigma = sim_sigma)
sim_graph <- ifelse(abs(sim_omega) > 0, 1, 0)
diag(sim_graph) <- 0
sum(sim_graph)
```

    ## [1] 1210

### 함수 적용 + 경과 시간 체크

``` r
# sssl
tic()
sssl_model <- concentration_sssl(Y, pi = 2/(p-1), burn = 500, iter = 500, h = 50)
```

    ## iter :  500 
    ## iter :  1000

``` r
toc()
```

    ## 113.28 sec elapsed

``` r
# horseshoe
tic()
horseshoe_model <- concentration_horseshoe(Y = Y, burn = 500, iter = 500, 
                                           approximate = FALSE)
```

    ## iter :  500 
    ## iter :  1000

``` r
toc()
```

    ## 120.93 sec elapsed

``` r
# approximate horseshoe with fixed threshold
tic()
fixed_approx_horseshoe_model <- concentration_horseshoe(Y = Y, burn = 500, 
                                                        iter = 500,
                                                        fixed_threshold = p^(-2),
                                                        approximate = TRUE)
```

    ## iter :  500 
    ## iter :  1000

``` r
toc()
```

    ## 54.87 sec elapsed

``` r
# approximate horseshoe with adaptive threshold
tic()
adaptive_approx_horseshoe_model <- concentration_horseshoe(Y = Y, burn = 500,
                                                           iter = 500,
                                                           approximate = TRUE, 
                                                           auto.threshold = TRUE)
```

    ## iter :  500 
    ## iter :  1000

``` r
toc()
```

    ## 58.53 sec elapsed

### Frobenius norm

``` r
# sssl
sum((sssl_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.8996065

``` r
# horseshoe
sum((horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.3246181

``` r
# approximate horseshoe with fixed threshold
sum((fixed_approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 1.024075

``` r
# approximate horseshoe with adaptive threshold
sum((adaptive_approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 1.134965

### MCC measure

``` r
graph_evaluation(sim_graph, sssl_model$Graph)
```

    ## MCC :  0.6601795 
    ## Accuracy :  0.9356 
    ## Sensitivity :  0.4677686 
    ## Specificity :  1 
    ## Precision :  1

``` r
graph_evaluation(sim_graph, fixed_approx_horseshoe_model$Graph)
```

    ## MCC :  0.2995612 
    ## Accuracy :  0.6476 
    ## Sensitivity :  0.8330579 
    ## Specificity :  0.6220705 
    ## Precision :  0.2327945

``` r
graph_evaluation(sim_graph, adaptive_approx_horseshoe_model$Graph)
```

    ## MCC :  0.2968319 
    ## Accuracy :  0.6454 
    ## Sensitivity :  0.831405 
    ## Specificity :  0.6197952 
    ## Precision :  0.2313707

# 2. sparsity 강도가 강한 데이터

precision matrix에 대한 시뮬레이션 데이터를 생성하고, 이에 대한 역행렬을
true precision matrix로 설정하였다.

``` r
N <- 10000
p <- 100
real_pi <- 0.01
sim_data <- simulate_sparse_matrix(real_pi, p)
sim_omega <- round(sim_data, digits = 3)
sim_sigma <- solve(sim_omega)
Y <- mvrnorm(n = N, mu = rep(0, p), Sigma = sim_sigma)
sim_graph <- ifelse(abs(sim_omega) > 0, 1, 0)
diag(sim_graph) <- 0
sum(sim_graph)
```

    ## [1] 132

### 함수 적용 + 경과 시간 체크

``` r
# sssl
tic()
sssl_model <- concentration_sssl(Y, pi = 2/(p-1), burn = 500, iter = 500, h = 50)
```

    ## iter :  500 
    ## iter :  1000

``` r
toc()
```

    ## 131.48 sec elapsed

``` r
# horseshoe
tic()
horseshoe_model <- concentration_horseshoe(Y = Y, burn = 500, iter = 500, 
                                           approximate = FALSE)
```

    ## iter :  500 
    ## iter :  1000

``` r
toc()
```

    ## 135.78 sec elapsed

``` r
# approximate horseshoe with fixed threshold
tic()
fixed_approx_horseshoe_model <- concentration_horseshoe(Y = Y, burn = 500, 
                                                        iter = 500,
                                                        fixed_threshold = p^(-2),
                                                        approximate = TRUE)
```

    ## iter :  500 
    ## iter :  1000

``` r
toc()
```

    ## 34.86 sec elapsed

``` r
# approximate horseshoe with adaptive threshold
tic()
adaptive_approx_horseshoe_model <- concentration_horseshoe(Y = Y, burn = 500,
                                                           iter = 500,
                                                           approximate = TRUE, 
                                                           auto.threshold = TRUE)
```

    ## iter :  500 
    ## iter :  1000

``` r
toc()
```

    ## 40.75 sec elapsed

### Frobenius norm

``` r
# sssl
sum((sssl_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.7040363

``` r
# horseshoe
sum((horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.05463534

``` r
# approximate horseshoe with fixed threshold
sum((fixed_approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.1261657

``` r
# approximate horseshoe with adaptive threshold
sum((adaptive_approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.07886219

### MCC measure

``` r
graph_evaluation(sim_graph, sssl_model$Graph)
```

    ## MCC :  0.7861809 
    ## Accuracy :  0.995 
    ## Sensitivity :  0.6212121 
    ## Specificity :  1 
    ## Precision :  1

``` r
graph_evaluation(sim_graph, fixed_approx_horseshoe_model$Graph)
```

    ## MCC :  0.5689095 
    ## Accuracy :  0.9806 
    ## Sensitivity :  0.8484848 
    ## Specificity :  0.9823672 
    ## Precision :  0.3916084

``` r
graph_evaluation(sim_graph, adaptive_approx_horseshoe_model$Graph)
```

    ## MCC :  0.7143626 
    ## Accuracy :  0.9908 
    ## Sensitivity :  0.8484848 
    ## Specificity :  0.9927037 
    ## Precision :  0.6086957
