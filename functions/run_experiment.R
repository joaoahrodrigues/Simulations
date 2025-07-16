#######################################
# Run experiment - Reference Element
#######################################
run_experiment <- function(n = 500, d, k, m, imbalance_ratio = 9, bins = NULL,
                  runs = 20, cont_method, cat_method, bin_method, farness = FALSE) {
  
  ########################################
  # Given a database, converts categorical variables into factors
  convert_to_factors <- function(df) {
  
  df <- df %>%
    mutate(across(starts_with("Cat_"), as.factor)) %>%
    mutate(Class = as.factor(Class))
  
  return(df)}
  
  ########################################
  # Yeo-Johnson function (used in farness)
  # Based on the internal function YJ from the cellWise package
  YJ <- function(y, lambda, stdToo = TRUE) {
    
    lowInd <- which(y < 0)
    highInd <- which(y >= 0)
    if (lambda != 0) {
      y[highInd] <- ((1 + y[highInd])^(lambda) - 1) / lambda}
    else {
      y[highInd] <- log(1 + y[highInd])}
    if (lambda != 2) {
      y[lowInd] <- -((1 - y[lowInd])^(2 - lambda) - 1) / (2 - lambda)}
    else {
      y[lowInd] <- -log(1 - y[lowInd])}
    if (stdToo) {
      if (length(y) > 1) {
        locScale <- cellWise::estLocScale(matrix(y, ncol = 1), type = "hubhub")
        zt <- (y - locScale$loc) / locScale$scale}
      else {
        zt <- y}}
    else {
      zt <- NULL}
    return(list(yt = y, zt = zt))}
  
  ########################################
  # Farness function (adapted - input distances to the reference element)
  Farness <- function(distMatrix) {
    
    indnz  <- which(distMatrix > 1e-10)
    farnz  <- distMatrix[indnz]
    farloc <- median(farnz, na.rm = TRUE)
    farsca <- mad(farnz, na.rm = TRUE)
    if (farsca < 1e-10) {farsca <- sd(farnz, na.rm = TRUE)}
    sfar <- scale(farnz, center = farloc, scale = farsca)
    YJout  <- cellWise::transfo(X = sfar, robust = TRUE, standardize = FALSE, checkPars = list(silent = TRUE))
    xt      <- YJout$Y
    tfarloc <- median(xt, na.rm = TRUE)
    tfarsca <- mad(xt, na.rm = TRUE)
    lambda  <- YJout$lambdahats
    origIndnz <- which(distMatrix > 1e-10)
    origFarnz <- distMatrix[origIndnz]
    zt        <- scale(YJ(scale(origFarnz, farloc, farsca), lambda = lambda, stdToo = FALSE)$yt, tfarloc, tfarsca)
    probs     <- rep(0, length(distMatrix))
    probs[origIndnz] <- pnorm(zt)
    
    return(farness = probs)}
  
  ########################################
  # Function to run each iteration
  run_once <- function(i) {
    
    # Generate dataset
    df <- generate_dataset(n, d, k, m, imbalance_ratio)
    df <- convert_to_factors(df)
    
    # Separate categorical and continuous part
    Xcat <- df[,1:k]
    Xcont <- df[,(k+1):(d+k)]
    
    numerical_vars <- colnames(df)[grepl("Cont_", colnames(df))]
    categorical_vars <- colnames(df)[grepl("Cat_", colnames(df))]
    categorical_vars <- c(categorical_vars, "Class")
    categorical_vars_no_class <- setdiff(categorical_vars, "Class")
    
    # Calculate reference element for numerical part
    if (cont_method == "mean") {
      num_values <- colMeans(Xcont)
    } else if (cont_method == "median") {
      num_values <- matrixStats::colMedians(as.matrix(Xcont))
    } else if (cont_method == "FastMCD") {
      fast_mcd_result <- CovMcd(Xcont)
      num_values <- fast_mcd_result@center
    } else {
      stop("Invalid cont_method")
    }
    
    # Calculate reference element for categorical part
    if (cat_method == "Smode") {
      reference_cat <- df %>% count(across(all_of(categorical_vars_no_class)), sort = TRUE) %>% slice(1)
    } else if (cat_method == "MCA") {
      mca_result <- MCA(Xcat, graph = FALSE)
      mca_coord <- mca_result$ind$coord # Coordinates in lower dimension
      centroid_mca <- colMedians(mca_coord) # Centroid in lower dimension
      distances <- apply(mca_coord, 1, function(row) sqrt(sum((row - centroid_mca)^2)))
      closest_idx <- which.min(distances) # Find the index of the observation that is the most similar to the centroid
      reference_cat <- df[closest_idx, categorical_vars_no_class] # Take the categorical part of that observation
    } else {
      stop("Invalid cat_method")
    }
    
    # Create the reference element 
    reference <- tibble(
      !!!setNames(reference_cat[1, 1:length(categorical_vars_no_class)], categorical_vars_no_class),
      !!!setNames(num_values, numerical_vars))

    for (col in categorical_vars_no_class) {
      reference[[col]] <- factor(reference[[col]], levels = levels(df[[col]]))
      df[[col]] <- factor(df[[col]], levels = levels(df[[col]]))}
    
    df_no_class <- df[, !names(df) %in% "Class"]
    df_with_ref <- rbind(df_no_class, reference) 
    
    dist_methods <- c("cond", "js", "gower", "huang")
    results <- list()
    
    for (method in dist_methods) {

      # The distances to the Reference Element appear in the last row (except the last entry)
      if (method %in% c("cond", "js")) {
        aux <- Ahmad.aux(x.cont = df_with_ref[, numerical_vars], x.cat = df_with_ref[, categorical_vars_no_class], 
              type = "Norm", bins = bins, only.categ = FALSE, js = (method == "js"), bin_method = bin_method)
        distances <- aux$dist[nrow(df_with_ref), -ncol(df_with_ref)] 
      } else {
        idcat <- 1:k
        idnum <- (k+1):(d+k)
        new_obs_dist <- as.matrix(distmix(df_with_ref, method = method, idnum = idnum, idcat = idcat))
        distances <- as.numeric(new_obs_dist[nrow(new_obs_dist), 1:(nrow(new_obs_dist) - 1)])}
      
      # If farness = TRUE, apply it to the distances 
      if (farness) distances <- Farness(distances)
      
      threshold <- quantile(distances, 0.9) 
      removed_indices <- which(distances >= threshold) 
      outlier_indices <- which(df$Class == 1) 
      acc <- length(intersect(removed_indices, outlier_indices)) / length(outlier_indices)
      
      outlier_labels <- rep(0, n)
      outlier_labels[outlier_indices] <- 1
      predicted_outliers <- ifelse(distances >= threshold, 1, 0) 
      
      TP <- sum(outlier_labels == 1 & predicted_outliers == 1) / sum(outlier_labels == 1)
      TN <- sum(outlier_labels == 0 & predicted_outliers == 0) / sum(outlier_labels == 0)
      FP <- sum(outlier_labels == 0 & predicted_outliers == 1) / sum(outlier_labels == 0)
      FN <- sum(outlier_labels == 1 & predicted_outliers == 0) / sum(outlier_labels == 1)
      
      recall_1 <- TP / (TP + FN)
      precision_0 <- TN / (TN + FN)
      recall_0 <- TN / (TN + FP)
      precision_1 <- TP / (TP + FP)
      
      epsilon <- 1e-8
      #f1_score <- ((2 * recall_1 * precision_1) / (recall_1 + precision_1 + epsilon) +
      #               (2 * recall_0 * precision_0) / (recall_0 + precision_0 + epsilon)) / 2
      f1_score <- (2 * recall_1 * precision_1) / (recall_1 + precision_1 + epsilon)
      balanced_accuracy <- (recall_0 + recall_1) / 2
      
      results[[method]] <- data.frame(
        Iteration = i, d = d, k = k, m = m, Method = method,
        Bins = if (method %in% c("cond", "js")) bins else NA,
        bin_method = if (method %in% c("cond", "js")) bin_method else NA,
        Accuracy = acc, F1_Score = f1_score, 
        Balanced_Accuracy = balanced_accuracy, Recall_1 = recall_1,
        Precision_0 = precision_0, cont_method = cont_method,
        cat_method = cat_method, farness = farness)}
    
    bind_rows(results)}

  # Run in parallel
  plan(multisession)  
  future_map_dfr(1:runs, run_once, .options = furrr_options(seed = TRUE))} # End function
