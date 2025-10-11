#######################################
# Generate dataset type 2
#######################################
generate_dataset_new <- function(n = 1000, d, k, m, imbalance_ratio = 9) {

  n_minority <- n / (1 + imbalance_ratio) # 10% of the data
  n_majority <- n - n_minority            # 90% of the data

  mu_majority <- rep(0, d)
  sigma_majority <- diag(d) * 1.5 # Variance of the dominant class

  generate_donut_values <- function(n, d, inner_radius = 3, outer_radius = 6) {

    result <- matrix(NA, nrow = n, ncol = d)
    for (i in 1:n) {
      repeat {
        candidate <- rnorm(d, mean = 0, sd = 1)
        r <- sqrt(sum(candidate^2))
        if (r >= inner_radius && r <= outer_radius) {
          result[i, ] <- candidate * (runif(1, inner_radius, outer_radius) / r)
          break}}}

    return(result)}

  X_con_minority <- generate_donut_values(n_minority, d, inner_radius = 3, outer_radius = 6)
  X_con_majority <- mvrnorm(n_majority, mu_majority, sigma_majority)

  generate_categorical <- function(n, k, m, class) {
    prob_matrix <- matrix(0, nrow = k, ncol = m)

    if (class == 0) {
      for (i in 1:k) {
        prob_matrix[i, 1] <- 0.4 + (i - 1) * (0.55 - 0.4) / (k - 1)
        prob_matrix[i, -1] <- rep((1 - prob_matrix[i, 1]) / (m - 1), m - 1)}
    } else {
      for (i in 1:k) {
        prob_matrix[i, m] <- 0.4 + (i - 1) * (0.55 - 0.4) / (k - 1)
        prob_matrix[i, 1:(m-1)] <- rep((1 - prob_matrix[i, m]) / (m - 1), m - 1)}}

    categorical_data <- data.frame(
      lapply(1:k, function(x) sample(1:m, n, replace = TRUE, prob = prob_matrix[x, ])))
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
