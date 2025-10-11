#######################################
# coocc_wt2 - compute the weights for the numerical variables
#######################################
coocc_wt2 <- function(x.cont, x.cat, x.disc, weights = NULL, bins = NULL, bin_method = NULL) {

  bin_method <- match.arg(bin_method, choices = c("equalwidth", "weighted"))
  if (is.null(weights)) weights <- rep(1, nrow(x.cont))

  wt <- rep(0, ncol(x.cont))

  for (h in 1:ncol(x.cont)) {
    df <- as.matrix(cbind(x.disc[, h], x.cat, x.disc[, -h, drop = FALSE]))
    delta <- 0

    for (j in 1:(ncol(df) - 1)) {
      tmp <- xtabs(~ df[, 1] + df[, -1][, j])
      tmp <- tmp / rowSums(tmp)
      delta <- delta + as.matrix(dist(tmp, "manhattan"))}

    delta <- 0.5 * delta / (ncol(df) - 1)
    wt[h] <- sum(delta) / ((ncol(delta) - 1) * ncol(delta))}

  return(wt)}
