Simulations
================
2024-01-08

### 만들어 둔 함수 불러오기

``` r
source(file = "./graph model/concentration_horseshoe_function.R")
source(file = "./graph model/concentration_sssl_function.R")
source(file = "./graph model/graph_evaluation_function.R")
source(file = "./graph model/simulate_data_function.R")
```

# concentration model

## 1. non-zero 10% 비율 시뮬레이션 데이터에 알고리즘 적용

p = 100으로 설정하고, concentration matrix에서 전체 가능한 edge 4950개
중 약 10%인 495개의 노드만이 활성화되어 있는 시뮬레이션 데이터를
생성하여, 다음 3가지 결과를 비교했다.

1.  norm 비교
2.  경과 시간 비교
3.  graph 결과 비교

### 시뮬레이션 데이터

``` r
N <- 10000
p <- 100
real_pi <- 0.1
sim_omega <- simulate_sparse_matrix(nonzero_ratio = real_pi, dimension = p)
sim_sigma <- solve(sim_omega)
Y <- mvrnorm(n = N, mu = rep(0, p), Sigma = sim_sigma)
```

### 함수 적용

``` r
tic()
sssl_model <- concentration_sssl(Y, pi = 2/(p-1), 
                                 burn = 500, iter = 1000, h = 50)
toc()
```

    ## 168.29 sec elapsed

``` r
tic()
horseshoe_model <- concentration_horseshoe(Y = Y, 
                                           burn = 500, iter = 1000, 
                                           approximate = FALSE)
toc()
```

    ## 179.49 sec elapsed

``` r
tic()
fix_approx_horseshoe_model <- concentration_horseshoe(Y = Y, 
                                                      burn = 500, iter = 1000,
                                                      fixed_threshold = 10^(-4),
                                                      approximate = TRUE)
toc()
```

    ## 47.78 sec elapsed

``` r
tic()
adapt_approx_horseshoe_model <- concentration_horseshoe(Y = Y, 
                                                        burn = 500, iter = 1000,
                                                        approximate = TRUE, 
                                                        auto.threshold = TRUE)
toc()
```

    ## 57.72 sec elapsed

### 추정된 precision matrix에 대한 확인

Omega 혹은 Sigma에 대한 true matrix를
$A=\{a_{ij} \}_{i,j=1}^{p} \in \mathbb{R}^{p \times p}$라고 표현할 때,

1.  sssl, horseshoe, approx_horseshoe 모델에 대하여 다음 norm을
    비교해보자.

$$||A - \hat{A}||_F^2$$

위의 값이 작을수록 좋은 추정을 했다고 판단할 수 있을 것이다.

2.  또한, 그려진 *그래프에 대해서 MCC 지표 등을 확인*하고자 한다.

### precision matrix에 대한 norm 확인

``` r
# sssl
sum((sssl_model$OmegaHat - sim_omega)^2)
```

    ## [1] 1.347521

``` r
# horseshoe
sum((horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.6840632

``` r
# fixed_approx_horseshoe
sum((fix_approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.4426328

``` r
# adapt_approx_horseshoe
sum((adapt_approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.7747205

### graph에 대한 지표 확인

``` r
true_graph <- ifelse(abs(sim_omega) > 0, 1, 0)
diag(true_graph) <- 0
sssl_graph <- sssl_model$Graph
fixed_approx_horseshoe_graph <- fix_approx_horseshoe_model$Graph
adapt_approx_horseshoe_graph <- adapt_approx_horseshoe_model$Graph

graph_evaluation(true_graph, sssl_graph)
```

    ## MCC :  1 
    ## Accuracy :  1 
    ## Sensitivity :  1 
    ## Specificity :  1 
    ## Precision :  1

``` r
graph_evaluation(true_graph, fixed_approx_horseshoe_graph)
```

    ## MCC :  1 
    ## Accuracy :  1 
    ## Sensitivity :  1 
    ## Specificity :  1 
    ## Precision :  1

``` r
graph_evaluation(true_graph, adapt_approx_horseshoe_graph)
```

    ## MCC :  0.9944152 
    ## Accuracy :  0.999 
    ## Sensitivity :  0.9979798 
    ## Specificity :  0.9991121 
    ## Precision :  0.9919679

## 2. non-zero 2% 비율 시뮬레이션 데이터에 알고리즘 적용

p = 100으로 설정하고, concentration matrix에서 전체 가능한 edge 4950개
중 약 2%인 99개의 노드만이 활성화되어 있는 시뮬레이션 데이터를 생성하여
동일한 방법으로 비교했다.

### 시뮬레이션 데이터 생성

``` r
N <- 10000
p <- 100
real_pi <- 0.02
sim_omega <- simulate_sparse_matrix(nonzero_ratio = real_pi, dimension = p)
sim_sigma <- solve(sim_omega)
Y <- mvrnorm(n = N, mu = rep(0, p), Sigma = sim_sigma)
```

### 함수 적용

``` r
tic()
sssl_model <- concentration_sssl(Y, pi = 2/(p-1), 
                                 burn = 500, iter = 1000, h = 50)
toc()
```

    ## 171.86 sec elapsed

``` r
tic()
horseshoe_model <- concentration_horseshoe(Y = Y, 
                                           burn = 500, iter = 1000, 
                                           approximate = FALSE)
toc()
```

    ## 175.82 sec elapsed

``` r
tic()
fix_approx_horseshoe_model <- concentration_horseshoe(Y = Y, 
                                                      burn = 500, iter = 1000,
                                                      fixed_threshold = 10^(-4),
                                                      approximate = TRUE)
toc()
```

    ## 42.43 sec elapsed

``` r
tic()
adapt_approx_horseshoe_model <- concentration_horseshoe(Y = Y, 
                                                        burn = 500, iter = 1000,
                                                        approximate = TRUE, 
                                                        auto.threshold = TRUE)
toc()
```

    ## 53.55 sec elapsed

### precision matrix에 대한 norm 확인

``` r
# sssl
sum((sssl_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.6464365

``` r
# horseshoe
sum((horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.03174351

``` r
# fixed_approx_horseshoe
sum((fix_approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.1848697

``` r
# adapt_approx_horseshoe
sum((adapt_approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 0.1447483

### graph에 대한 지표 확인

``` r
true_graph <- ifelse(abs(sim_omega) > 0, 1, 0)
diag(true_graph) <- 0
sssl_graph <- sssl_model$Graph
fixed_approx_horseshoe_graph <- fix_approx_horseshoe_model$Graph
adapt_approx_horseshoe_graph <- adapt_approx_horseshoe_model$Graph

graph_evaluation(true_graph, sssl_graph)
```

    ## MCC :  1 
    ## Accuracy :  1 
    ## Sensitivity :  1 
    ## Specificity :  1 
    ## Precision :  1

``` r
graph_evaluation(true_graph, fixed_approx_horseshoe_graph)
```

    ## MCC :  1 
    ## Accuracy :  1 
    ## Sensitivity :  1 
    ## Specificity :  1 
    ## Precision :  1

``` r
graph_evaluation(true_graph, adapt_approx_horseshoe_graph)
```

    ## MCC :  0.9611053 
    ## Accuracy :  0.9984 
    ## Sensitivity :  1 
    ## Specificity :  0.9983677 
    ## Precision :  0.9252336

## 3. non-zero 50% 비율 시뮬레이션 데이터에 알고리즘 적용

해당 시뮬레이션은 non-zero 비율이 높을 때의 결과를 확인하기 위함이다.

### 시뮬레이션 데이터 생성

``` r
N <- 10000
p <- 100
real_pi <- 0.5
sim_omega <- simulate_sparse_matrix(nonzero_ratio = real_pi, dimension = p)
sim_sigma <- solve(sim_omega)
Y <- mvrnorm(n = N, mu = rep(0, p), Sigma = sim_sigma)
```

### 함수 적용

``` r
tic()
sssl_model <- concentration_sssl(Y, pi = 2/(p-1), 
                                 burn = 500, iter = 1000, h = 50)
toc()
```

    ## 167.05 sec elapsed

``` r
tic()
horseshoe_model <- concentration_horseshoe(Y = Y, 
                                           burn = 500, iter = 1000, 
                                           approximate = FALSE)
toc()
```

    ## 174.36 sec elapsed

``` r
tic()
fix_approx_horseshoe_model <- concentration_horseshoe(Y = Y, 
                                                      burn = 500, iter = 1000,
                                                      fixed_threshold = 10^(-4),
                                                      approximate = TRUE)
toc()
```

    ## 146.65 sec elapsed

``` r
tic()
adapt_approx_horseshoe_model <- concentration_horseshoe(Y = Y, 
                                                        burn = 500, iter = 1000,
                                                        approximate = TRUE, 
                                                        auto.threshold = TRUE)
toc()
```

    ## 134.08 sec elapsed

### precision matrix에 대한 norm 확인

``` r
# sssl
sum((sssl_model$OmegaHat - sim_omega)^2)
```

    ## [1] 11.0101

``` r
# horseshoe
sum((horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 18.37482

``` r
# fixed_approx_horseshoe
sum((fix_approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 16.75782

``` r
# adapt_approx_horseshoe
sum((adapt_approx_horseshoe_model$OmegaHat - sim_omega)^2)
```

    ## [1] 14.68189

### graph에 대한 지표 확인

``` r
true_graph <- ifelse(abs(sim_omega) > 0, 1, 0)
diag(true_graph) <- 0
sssl_graph <- sssl_model$Graph
fixed_approx_horseshoe_graph <- fix_approx_horseshoe_model$Graph
adapt_approx_horseshoe_graph <- adapt_approx_horseshoe_model$Graph

graph_evaluation(true_graph, sssl_graph)
```

    ## MCC :  1 
    ## Accuracy :  1 
    ## Sensitivity :  1 
    ## Specificity :  1 
    ## Precision :  1

``` r
graph_evaluation(true_graph, fixed_approx_horseshoe_graph)
```

    ## MCC :  0.1229999 
    ## Accuracy :  0.5102 
    ## Sensitivity :  1 
    ## Specificity :  0.03009901 
    ## Precision :  0.5026401

``` r
graph_evaluation(true_graph, adapt_approx_horseshoe_graph)
```

    ## MCC :  0.486967 
    ## Accuracy :  0.6898 
    ## Sensitivity :  1 
    ## Specificity :  0.3857426 
    ## Precision :  0.6147541
