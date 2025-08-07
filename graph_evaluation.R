graph_evaluation <- function(true_graph, pred_graph) {
  TP <- sum((true_graph == 1) & (pred_graph == 1))
  FP <- sum((true_graph == 0) & (pred_graph == 1))
  TN <- sum((true_graph == 0) & (pred_graph == 0))
  FN <- sum((true_graph == 1) & (pred_graph == 0))
  
  MCC <- (TP * TN - FP * FN) / sqrt(TP + FP)
  MCC <- MCC / sqrt(TP + FN)
  MCC <- MCC / sqrt(TN + FP)
  MCC <- MCC / sqrt(TN + FN)
  
  Accuracy <- (TP + TN)/(TP + TN + FP + FN)
  Sensitivity <- TP/(TP + FN)
  Specificity <- TN/(TN + FP)
  Precision <- TP/(TP + FP)
  
  cat("MCC : ", MCC, "\n")
  cat("Accuracy : ", Accuracy, "\n")
  cat("Sensitivity : ", Sensitivity, "\n")
  cat("Specificity : ", Specificity, "\n")
  cat("Precision : ", Precision, "\n")
}