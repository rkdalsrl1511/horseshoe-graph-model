# 시뮬레이션 데이터 생성 함수
simulate_sparse_matrix <- function(nonzero_ratio, dimension) {
  
  if (abs(nonzero_ratio) > 1 | dimension < 5) {
    stop("0<=nonzero_ratio<=1, dimension >= 5를 입력하라.")
  }
  result_matrix <- diag(1, nrow = dimension, ncol = dimension)
  n_nonzero <- ceiling(nonzero_ratio * dimension * (dimension-1) / 2)
  graph_list <- list()
  k <- 1
  for (i in 2:dimension) {
    for (j in 1:(i-1)) {
      graph_list[[k]] <- c(j, i)
      k <- k + 1
    }
  }
  # non-zero element's indexes
  index_nonzero <- sample(graph_list, size = n_nonzero)
  # Random sampling of non-zero element values
  for (i in 1:length(index_nonzero)) {
    temp <- index_nonzero[[i]]
    v <- sample(x = c(-1, 1), size = 1)
    u <- runif(1, 0.3, 0.5)
    result_matrix[temp[1], temp[2]] <- v * u
    result_matrix[temp[2], temp[1]] <- v * u
  }
  # positive definite 맞추기
  min_eigen_value <- min(eigen(result_matrix)$value)
  if (min_eigen_value <= 0) {
    eigen_diag <- diag(abs(min_eigen_value) + 0.001, 
                       nrow = dimension, ncol = dimension)
    result_matrix <- result_matrix + eigen_diag
  }
  round(result_matrix, digits = 3)
}