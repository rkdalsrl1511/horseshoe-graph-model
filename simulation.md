simulation
================
2024-04-09

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.2 ──
    ## ✔ ggplot2 3.4.1      ✔ purrr   0.3.4 
    ## ✔ tibble  3.2.1      ✔ dplyr   1.0.10
    ## ✔ tidyr   1.2.0      ✔ stringr 1.4.1 
    ## ✔ readr   2.1.2      ✔ forcats 0.5.2

    ## Warning: 패키지 'ggplot2'는 R 버전 4.2.3에서 작성되었습니다

    ## Warning: 패키지 'tibble'는 R 버전 4.2.3에서 작성되었습니다

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
library(MASS)
```

    ## 
    ## 다음의 패키지를 부착합니다: 'MASS'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

``` r
library(tictoc)
library(GIGrvg)
```

    ## Warning: 패키지 'GIGrvg'는 R 버전 4.2.3에서 작성되었습니다

### 시뮬레이션 데이터 생성

``` r
N <- 10000
p <- 200
real_pi <- 0.05
sim_data <- simulate_sparse_matrix(real_pi, p)
sim_sigma <- round(sim_data, digits = 3)
sim_omega <- solve(sim_sigma)
Y <- mvrnorm(n = N, mu = rep(0, p), Sigma = sim_sigma)
```

### 함수 적용 + 경과 시간 체크

``` r
tic() 
approx_horseshoe_model <- concentration_horseshoe(Y = Y, burn = 0, iter = 500,
                                                  fixed_threshold = 10^(-4),
                                                  approximate = TRUE)
```

    ## iter :  500

``` r
toc()
```

    ## 266.63 sec elapsed

``` r
tic()
auto_approx_horseshoe_model <- concentration_horseshoe(Y = Y, burn = 0, 
                                                       iter = 500,
                                                       approximate = TRUE, 
                                                       auto.threshold = TRUE)
```

    ## iter :  500

``` r
toc()
```

    ## 279.3 sec elapsed

``` r
tic()
sssl_model <- concentration_sssl(Y, pi = 2/(p-1), burn = 0, iter = 500, h = 50)
```

    ## iter :  500

``` r
toc()
```

    ## 512.93 sec elapsed

``` r
# sssl
sum((sssl_model$OmegaHat - sim_omega)^2)
```

    ## [1] 95.28329

``` r
# approx_horseshoe
sum((approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 48.84791

``` r
# auto_approx_horseshoe
sum((auto_approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 44.67175
