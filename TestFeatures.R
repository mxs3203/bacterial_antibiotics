#!/usr/bin/env Rscript

#### Load Needed packages, environment and save parameters ####
suppressPackageStartupMessages({
    library(tidyverse)
    library(tidyr)
})

cat("\nLoading .RData file...\n")

load("ARG.RData")

args = commandArgs(trailingOnly=TRUE)

iterations <- args[1]

#### TEST FOR SIGNIFICANT COLUMNS FUNCTION ####
test_cols <- function(all.genes, iterations){
    
    resistent.i <- which(all.genes["resistance"] == "T")    # take row number of resistance genes
    non_resistent.i <- which(all.genes["resistance"] == "F")
    
    results.list <- lapply(all.genes[,12:354], function(x){ # pass each column at a time
        totalMedian = 0
        totalMean = 0
        t <- replicate(iterations, {    
            # Resample resistent and non-resistent genes data from passed column
            resampled_resistent <- sample(x[resistent.i], length(resistent.i), replace = T)
            resampled_non_resistent <- sample(x[non_resistent.i], 1000, replace = T)
            
            test_results <- t.test(resampled_resistent, resampled_non_resistent)
            
            totalMedian = totalMedian + median(resampled_resistent) - median(resampled_non_resistent)
            totalMean = totalMean + mean(resampled_resistent) - mean(resampled_non_resistent)
            
            return(c("totalMean" = totalMean, "totalMedian" = totalMedian,
                     test_results$statistic, p.value = test_results$p.value))   # Returns results in a new column of t
        })
        
        t <- apply(t(t), 2, function(x) abs(mean(x)))   # Calculate the abs mean of the results from n iterations
        
    })
    
    return(do.call("rbind", results.list))  # bind every list element in a data.frame with element names as row names
}

cat("\nTesting feature differences between resistant and non-resistance genes by two sample t - test
    Nº iterations: ", iterations, "\n\n")

#### TEST USING FUNCTION ####
significant_cols_results <- test_cols(all.genes, iterations)
significant_cols_results <- as.data.frame(significant_cols_results)
significant_cols_results <- mutate(significant_cols_results, col = rownames(significant_cols_results))

cat("First rows of the results, arranged by p.value:\n")
print(significant_cols_results %>%
    arrange(p.value) %>%
    head())

cat("\nSaving results to 'ARG.RData'...\n")

save.image(file = "ARG.RData")

cat("\n\t'ARG.RData' file updated :)\n\n")

