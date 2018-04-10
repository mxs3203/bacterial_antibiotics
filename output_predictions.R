args = commandArgs(trailingOnly=TRUE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(tidyr)
  library(stats)
  library(glmnet)
  library(class)
})

method <- args[1]
fileName <- args[2]


#method <- "GLM"
if(method == "KNN") {
  load(file='Models/ARG_KNN.RData')
} else if(method == "GLM"){
  load(file='Models/ARG_GLM.RData')
}


new_data = read.csv(fileName, header = T)

ten_signficant_vars <- significant_cols_results %>%
    arrange(p.value) %>%
    head(number_explanatory_var)

new_data <- new_data[ , which(names(new_data) %in% c("resistance", ten_signficant_vars$col))] 

test <- model.matrix(resistance~., new_data, type = "response")
response <- predict(trained_model, newx = test, s = "lambda.min", type = "class")
write.csv(response, file = "predicted.csv")

cat("\nPredictions are written in file: predicted.csv \n \n ")
