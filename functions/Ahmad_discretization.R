#######################################
# generate_binned_cont -> discretize numerical variables
#######################################
generate_binned_cont <- function(x.cont, bins = NULL, bin_method = NULL) {
  
  bin_method <- match.arg(bin_method, choices = c("equalwidth", "weighted"))
  if (is.null(bin_method)) bin_method <- "equalwidth"
  if (is.null(bins)) bins <- ceiling(sqrt(nrow(x.cont)))
  # Compute Mahalanobis-based weights
  if (bin_method == "weighted") {
    mcd_fit <- covMcd(x.cont)
    mahal_dist <- mahalanobis(x.cont, center = mcd_fit$center, cov = mcd_fit$cov)
    weights <- 1 / (1 + mahal_dist)}
  
  weighted_bin_disc <- function(x, weights, bins) {
    probs <- seq(0, 1, length.out = bins + 1)
    cut_points <- wtd.quantile(x, weights = weights, probs = probs, na.rm = TRUE)
    cut_points <- sort(unique(cut_points))
    
    if (length(cut_points) < 2 || any(is.na(cut_points))) {
      cut_points <- seq(min(x, na.rm = TRUE) - 1e-8, max(x, na.rm = TRUE) + 1e-8, length.out = bins + 1)
    } else {
      eps <- .Machine$double.eps * 100
      cut_points[1] <- cut_points[1] - eps
      cut_points[length(cut_points)] <- cut_points[length(cut_points)] + eps}
    
    bins <- cut(x, breaks = cut_points, include.lowest = TRUE, labels = FALSE, right = FALSE)
    bins[is.na(bins) & x <= min(cut_points)] <- 1
    bins[is.na(bins) & x >= max(cut_points)] <- length(cut_points) - 1
    return(as.factor(bins))}
  
  binned_cont <- as.data.frame(lapply(x.cont, function(col) {
    if (bin_method == "weighted") {
      weighted_bin_disc(x = col, weights = weights, bins = bins)
    } else {
      cut(col, breaks = bins, include.lowest = TRUE, labels = FALSE) %>% as.factor()}}))
  
  return(binned_cont)}
