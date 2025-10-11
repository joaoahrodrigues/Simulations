#######################################
# Ahmad.aux - adapted with Mahalanobis
#######################################
AhmadMah.aux <- function(x.cont, x.cat, type = "Norm", bins = NULL, only.categ = FALSE, js = FALSE, bin_method = NULL, method_mah = NULL) {

  method_mah <- match.arg(method_mah, choices = c("st", "rb"))
  bin_method <- match.arg(bin_method, choices = c("equalwidth", "weighted"))
  x.disc <- generate_binned_cont(as.data.frame(x.cont), bins = bins, bin_method = bin_method)
  d.cat <- coocc_dis(x.cat, x.cont, x.disc, bins = bins, only.cat = only.categ, js = js, bin_method = bin_method)

  if (type == "MinMax") {
    process <- caret::preProcess(as.data.frame(x.cont), method = c("range"))
    x.conts <- predict(process, as.data.frame(x.cont))
    wt <- coocc_wt2(x.conts, x.cat, x.disc, bins = bins, bin_method = bin_method)
    newdf <- sweep(x.conts, 2, wt^2, function(x, y) x * sqrt(y))

  } else if (type == "Norm") {
    wt <- coocc_wt2(scale(x.cont), x.cat, x.disc, bins = bins, bin_method = bin_method)
    newdf <- sweep(scale(x.cont), 2, wt^2, function(x, y) x * sqrt(y))

  } else if (type == "None") {
    wt <- rep(1, ncol(x.cont))
    newdf <- sweep(x.cont, 2, wt^2, function(x, y) x * sqrt(y))}

  mahal_dist_matrix <- function(X, center, cov) {

    Xc <- sweep(X, 2, center)
    inv_cov <- solve(cov)
    Q <- Xc %*% inv_cov %*% t(Xc)
    d2 <- outer(diag(Q), diag(Q), "+") - 2*Q
    d2[d2 < 0] <- 0

    return(sqrt(d2))}   

  newdf <- as.matrix(newdf)

  if (method_mah == "rb") {
    mcd_fit <- covMcd(newdf)
    d.cont <- mahal_dist_matrix(newdf, mcd_fit$center, mcd_fit$cov)^2

  } else { # standard method
    center <- colMeans(newdf)
    cov_mat <- cov(newdf)
    d.cont <- mahal_dist_matrix(newdf, center, cov_mat)^2}

  d.Ahmad <- d.cat + d.cont

  return(list(dist=sqrt(d.Ahmad), cont=d.cont, cat=d.cat))}
