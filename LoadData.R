#!/usr/bin/env Rscript

#### Load Needed packages and save parameters ####

suppressPackageStartupMessages({
    library(tidyverse)
    library(tidyr)
})

args = commandArgs(trailingOnly=TRUE)
csv_dir <- args[1]
# csv_dir <- "Data-Csv"
cat("CSV directory: ", path.expand(csv_dir), "\n\n")

#### Load CSV files ####

# Save csv file name as data frame name
temp <- list.files(path = csv_dir, pattern="*.csv")
df_names <- gsub("standardized", "", gsub("*.csv$", "", temp)) # remove "standardized" and .csv extension"
df_names <- make.names(df_names)    # replace problematic characters

cat("Reading CSV files...   ")

# Load each csv file, set the proper name for it and store each in a variable (data.frame)
temp <- list.files(path = csv_dir, pattern="*.csv", full.names = T)
suppressWarnings(suppressMessages(
    list2env(lapply(setNames(temp, df_names), read_csv), envir = .GlobalEnv)
))

cat("\n\tCSV files Loaded\n")
cat("\nPreprocessing...\n")


#### Preprocessing: Filters and get strain+specie columns ####
Ecoli_IAI39 <- Ecoli_IAI39 %>%
    mutate(strain = "IAI39", resistance = "F")
Ecoli_K12 <- Ecoli_K12 %>%
    mutate(strain = "K12", resistance = "F")
EColi_O104.H4 <- EColi_O104.H4 %>%
    mutate(strain = "O104.H4", resistance = "F")
EColi_O157.H7 <- EColi_O157.H7 %>%
    mutate(strain = "O157.H7", resistance = "F")
EColi_UMN026 <- EColi_UMN026 %>%
    mutate(strain = "UMN026", resistance = "F")
EColi_O83.H1 <- EColi_O83.H1 %>%
    mutate(strain = "O83.H1", resistance = "F")

ShigellaDysenteriae <- ShigellaDysenteriae %>%
    mutate(strain = "Dysenteriae", resistance = "F")
ShigellaDysenteriae$specie <- "Shigella"

ShigellaFlexneri <- ShigellaFlexneri %>%
    mutate(strain = "Flexneri", resistance = "F")
ShigellaFlexneri$specie <- "Shigella"

test <- finalResistantGenes

test$strain <- if_else(grepl("Escherichia coli", x = test$specie), gsub("Escherichia coli", "", test$specie), 
                       test$specie)
test$strain <- if_else(grepl("Shigella", x = test$specie), gsub("Shigella", "", test$specie), 
                       test$strain)

test$strain <- if_else(grepl("K-12", x = test$strain), "K12", 
                       if_else(grepl("flexneri", x = test$strain), "Flexneri", 
                               if_else(grepl("O157:H7", x = test$strain), "O157.H7", 
                                       gsub(pattern = "^ ", replacement = "", x = test$strain))))

test$specie <- if_else(grepl("Escherichia coli", x = test$specie), "Escherichia coli", 
                       if_else(grepl("Shigella", x = test$specie), "Shigella", test$specie))

finalResistantGenes <- test %>% mutate(resistance = "T")

all.genes <- rbind(Ecoli_IAI39, Ecoli_K12,EColi_O104.H4,EColi_O157.H7, 
                   EColi_UMN026, EColi_O83.H1, ShigellaDysenteriae,ShigellaFlexneri, finalResistantGenes)

cat("\nSaving Results to 'ARG.RData'...\n")

save.image(file = "ARG.RData")

cat("\n\t'ARG.RData' file saved :)\n\n")
