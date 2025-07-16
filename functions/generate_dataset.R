#######################################
# Generate dataset function 
#######################################
generate_dataset <- function(n = 1000, d, k, m, imbalance_ratio = 9) {
  
  n_minority <- n / (1 + imbalance_ratio) # 100 observations in minority class
  n_majority <- n - n_minority # 900 observations in majority class
  mu_majority <- rep(1.5, d) # Mean majority class
  mu_minority <- rep(-1.5, d) # Mean minority class
  sigma <- diag(d) * 1.5 # Diagonal covariance matrix with variance 1.5
  
  # Generate continuous variables (Normal distribution)
  X_con_majority <- mvrnorm(n_majority, mu_majority, sigma) 
  X_con_minority <- mvrnorm(n_minority, mu_minority, sigma)
  
  # Function to generate categorical variables
  generate_categorical <- function(n, k, m, class) {
    
    prob_matrix <- matrix(0, nrow = k, ncol = m)
    if (class == 0) {
      for (i in 1:k) {
        prob_matrix[i, 1] <- 0.4 + (i - 1) * (0.55 - 0.4) / (k - 1)
        prob_matrix[i, -1] <- rep((1 - prob_matrix[i, 1]) / (m - 1), m - 1)
        cat("Class 0 Probability Matrix - Variable", i, "\n", prob_matrix[i, ], "\n")}
    } else {
      for (i in 1:k) {
        prob_matrix[i, m] <- 0.4 + (i - 1) * (0.55 - 0.4) / (k - 1)
        prob_matrix[i, 1:(m-1)] <- rep((1 - prob_matrix[i, m]) / (m - 1), m - 1)
        cat("Class 1 Probability Matrix - Variable", i, "\n", prob_matrix[i, ], "\n")}}
    categorical_data <- data.frame(
      lapply(1:k, function(x) sample(1:m, n, replace = TRUE, prob = prob_matrix[x,])))
    colnames(categorical_data) <- paste("Cat", 1:k, sep = "_")
    
    return(categorical_data)}
  
  X_cat_majority <- generate_categorical(n_majority, k, m, class = 0)
  X_cat_minority <- generate_categorical(n_minority, k, m, class = 1)
  df_majority <- data.frame(X_cat_majority, as.data.frame(X_con_majority), Class = 0)
  df_minority <- data.frame(X_cat_minority, as.data.frame(X_con_minority), Class = 1)
  colnames(df_majority)[(k+1):(k+d)] <- paste("Cont", 1:d, sep = "_")
  colnames(df_minority)[(k+1):(k+d)] <- paste("Cont", 1:d, sep = "_")
  df <- rbind(df_majority, df_minority)
  
  return(df)}
