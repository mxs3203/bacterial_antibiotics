#!/usr/bin/env Rscript

# setwd("~/Desktop/AarhusUniversity/resistencetoantibiotics")
# use this to install everything and then use that environment in your IDE
#conda install -c r rpy2

# Load libraries
suppressPackageStartupMessages({
    library(tidyverse)
    library(caret)
    library(MASS)
    library(purrr)
    library(tidyr)
    library(earth)
    library(ggplot2)
    library(glmnet)
    library(class) #knn
})


# Save arguments
args = commandArgs(trailingOnly=TRUE)

method = toupper(as.character(args[1]))
iterationsForModel = as.numeric(args[2])
number_explanatory_var = as.numeric(args[3])
KNN_Tune_len = as.numeric(args[4])
KNN_CV_repeats = as.numeric(args[5])

cat("\nLoading .RData... ")

load(file='ARG.RData')

cat("\n\nBuilding model with following parameters:\n")
cat("\n\tMethod -> ", method)
cat("\n\tIterations For Model -> ", iterationsForModel)
cat("\n\tNº of significant parameters used -> ", number_explanatory_var)
if (method == "KNN") {
    cat("\n\tKNN tune length -> ", KNN_Tune_len)
    cat("\n\tKNN Cross Validations repeats -> ", KNN_CV_repeats)
}else {
    if (method == "GLM") {
        cat("\n\nIgnoring parameters: KNN_Tune_len, KNN_CV_repeats" )
    }
}
cat("\n\n")


makeModel <- function(iterationsForModel = 40, iterationForSignificantColumns = 100, method = "GLM", parametersUsedForModel = -1, explanatoryVar_Pval_cutoff = 0.001,number_explanatory_var = 10, data = all.genes, responseVariable = "resistance", KNN_Tune_len = 15, KNN_CV_repeats = 5){
    
    # The significant_cols_results should be already in environment
    if (parametersUsedForModel == -1){
        cat("Significant paramater names (columns) not specifed.\nUsing ", number_explanatory_var, " most significant ones from two sample t - test\n\n")
        ten_signficant_vars <- significant_cols_results %>%
            # filter(p.value < explanatoryVar_Pval_cutoff) %>%
            arrange(p.value) %>%
            head(number_explanatory_var)
    } else {
        cat("Significant Paramaters (columns) specifed: Ignoring parameters explanatory, Var_Pval_cutoff, number_explanatory_var\n\n")
        ten_signficant_vars <- parametersUsedForModel
    }
    
    #DataFrame for return in order to visualize results
    results <- data.frame(falseNegatives=integer(),falsePositives=integer(),total=integer(), percantage = double(), Sensitivty=double(), Precision=double(), BalancedAccuracy=double(), AccuracyPVal=double, McnemarPValue=double())
    
    cat("Using method -> " , method, ": Starting ", iterationsForModel, " iterations \n\n")
    cat("Significant Paramaters: \n", ten_signficant_vars$col, "\n\n")
    
    for(i in 1:iterationsForModel){
        
        resistent <- data %>%
            filter(resistance == "T")
        non_resistent <- data %>%
            filter(resistance == "F")
        
        # model_resistent <- sample_n(filter(resistent, !is.na(resistance)), 390, replace = T)
        # model_non_resistent <- sample_n(filter(non_resistent, !is.na(resistance)), 2000, replace = T)
        # 
        # train_resistent_ind <- sample(seq_len(nrow(model_resistent)), 
        #                               size = floor(0.80 * nrow(model_resistent)))
        # train_non_resistent_ind <- sample(seq_len(nrow(model_non_resistent)), 
        #                                   size = floor(0.80 * nrow(model_non_resistent)))
        # 
        # train <- rbind(model_resistent[train_resistent_ind, ], model_non_resistent[train_non_resistent_ind, ])
        # test <- rbind(model_resistent[-train_resistent_ind, ], model_non_resistent[-train_non_resistent_ind, ])
        
        
        if (method == "KNN"){
            train_resistent <- sample_n(filter(resistent, !is.na(resistance)), 398, replace = T)
            train_resistent_80 <- sample_n(train_resistent, nrow(train_resistent)*0.8)
            train_non_resistent <- sample_n(filter(non_resistent, !is.na(resistance)), 1800, replace = T)
            train_non_resistent_80 <- sample_n(train_non_resistent, nrow(train_non_resistent)*0.8)
            
            tmp_test <- rbind(train_non_resistent, train_resistent)
            train <- rbind(train_non_resistent_80, train_resistent_80)
            
            test <- suppressMessages(
                anti_join(tmp_test, train)
            )
            
            
            #Use only specified columns and response variable
            train <- train[ , which(names(train) %in% c(responseVariable, ten_signficant_vars$col))] 
            test <- test[ , which(names(test) %in% c(responseVariable, ten_signficant_vars$col))] 
            
            #remove NAs if there is any...
            train <- train[,colSums(is.na(train))<nrow(train)]
            train <- na.omit(train)
            
            ctrl <- trainControl(method="repeatedcv", repeats = KNN_CV_repeats)
            trained_model <- train(resistance ~ . , data = train, method = "knn", trControl = ctrl, preProcess = c("center","scale"), tuneLength = KNN_Tune_len)
            
            
            response <- predict(trained_model, newdata = test)
            
            
        } else if (method == "GLM"){
            train_resistent <- sample_n(filter(resistent, !is.na(resistance)), 398, replace = T)
            train_resistent_80 <- sample_n(train_resistent, nrow(train_resistent)*0.8)
            train_non_resistent <- sample_n(filter(non_resistent, !is.na(resistance)), 1050, replace = T)
            train_non_resistent_80 <- sample_n(train_non_resistent, nrow(train_non_resistent)*0.8)
            
            tmp_test <- rbind(train_non_resistent, train_resistent)
            train <- rbind(train_non_resistent_80, train_resistent_80)
            
            test <- suppressMessages(
                anti_join(tmp_test, train)
            )
            
            
            #Use only specified columns and response variable
            train <- train[ , which(names(train) %in% c(responseVariable, ten_signficant_vars$col))] 
            test <- test[ , which(names(test) %in% c(responseVariable, ten_signficant_vars$col))] 
            
            #remove NAs if there is any...
            train <- train[,colSums(is.na(train))<nrow(train)]
            train <- na.omit(train)
            
            
            x <- model.matrix(resistance ~ . , data = train)
            y <- train$resistance
            #Cross validating lambda min
            trained_model <- cv.glmnet(x, y, family = "binomial", type.measure = "class")
            test1 <- model.matrix(resistance ~ ., train)
            #Remove NAs if there is any
            test <- test[,colSums(is.na(test))<nrow(test)]
            test <- na.omit(test)
            
            test2 <- model.matrix(resistance~., test)
            response <- predict(trained_model, newx = test2,
                                s = "lambda.min",
                                type = "class")

            
        } else {
            cat("Unknown method please use KNN or GLM \n\n")
            return -1
        }
        
        
        trained_model <<- trained_model
        
        confusion_matrix <- table(response, test$resistance)
        false_negatives <- confusion_matrix[1,2]
        false_positives <- confusion_matrix[2,1]
        total_number_of_wrong <- false_negatives + false_positives
        
        cm <- confusionMatrix(confusion_matrix)
        
        t = data.frame(false_negatives,false_positives,total_number_of_wrong, total_number_of_wrong/nrow(test),cm$byClass['Sensitivity'], cm$byClass['Precision'], cm$byClass['Balanced Accuracy'],cm$overall['AccuracyPValue'], cm$overall['McnemarPValue']) 
        names(t)=c("falseNegatives","falsePositives","total","percentage", "Sensitivty", "Precision", "BalancedAccuracy", "AccuracyPVal", "McnemarPValue") 
        results = rbind(results, t)
        
    }
    return (results)
}


results_df <- makeModel(iterationsForModel = iterationsForModel, method = method, 
                        number_explanatory_var = number_explanatory_var,
                        KNN_Tune_len = KNN_Tune_len, KNN_CV_repeats = KNN_CV_repeats)
cat("\nSummary of the results from testing the model:\n")
print(summary(results_df))
hist(results_df$AccuracyPVal)
ggplot(data = results_df) +
    geom_histogram(mapping = aes(x = AccuracyPVal), bins = 8, binwidth = 0.5)
cat("\nSaving results to 'ARG.RData'...\n")

suppressWarnings(rm(iterations, iterationsForModel, method, KNN_Tune_len, KNN_CV_repeats, args))

save.image(file = "ARG.RData")

cat("\n\t'ARG.RData' file updated :)\n\n")
