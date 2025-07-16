###############################################
# Master script to run the experiment
###############################################

# Load required libraries
library(robustbase)
library(cluster)  
library(ggplot2)  
library(factoextra)  
library(clValid)  
library(fpc)  
library(gower)
library("tidyverse")
library("haven")
library("dplyr")
library("cvTools")
library(boot)
library(textmineR)
library(caret)
library(proxy)
library(Matrix) 
library(expm) 
library(FNN) 
library(Rlof) 
library(dbscan) 
library(kmed)
library(data.table)
library(dplyr)
library(stringr)
library(parallel)
library(gridExtra)
library(pROC)
library(reshape2)
library(mvtnorm)
library(MASS)
library(purrr)
library(knitr)  
library(gt)
library(Rcpp)
library(rrcov)
library(FactoMineR)
library(factoextra)
library(philentropy)
library(furrr)
library(tibble)
library(Hmisc)
library(klaR)
library(progressr)
library(future)

# Source all functions from the 'functions' directory
function_files <- list.files("functions", full.names = TRUE, pattern = "\\.R$")
sapply(function_files, source)

# Set parameters
set.seed(123)
n <- 1000
d <- 3
k <- 3
m <- 3
imbalance_ratio <- 9
bins <- 20

# Example: Generate dataset and compute distances
df <- generate_dataset(n = n, d = d, k = k, m = m, imbalance_ratio = imbalance_ratio)
x.cat <- df[, 1:k]
x.cont <- df[, (k+1):(k+d)]

# Run Ahmad auxiliary distance
result <- Ahmad.aux(
  x.cont = x.cont,
  x.cat = x.cat,
  type = "MinMax",
  bins = bins,
  only.categ = FALSE,
  js = TRUE,
  bin_method = "equalwidth"
)

distances <- result$dist
print(distances[1:5, 1:5])  # Example output
