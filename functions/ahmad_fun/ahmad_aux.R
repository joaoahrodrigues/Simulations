#######################################
# Ahmad.aux
#######################################
Ahmad.aux <- function(x.cont, x.cat, type = "Norm", bins = NULL, only.categ = FALSE, js = FALSE, bin_method = NULL) {

  bin_method <- match.arg(bin_method, choices = c("equalwidth", "weighted"))
  x.disc <- generate_binned_cont(as.data.frame(x.cont), bins = bins, bin_method = bin_method)
  d.cat <- coocc_dis(x.cat, x.cont, x.disc, bins = bins, only.cat = only.categ, js = js, bin_method = bin_method)

  if (type == "MinMax") {
    process <- caret::preProcess(as.data.frame(x.cont), method = c("range"))
    x.conts <- predict(process, as.data.frame(x.cont))
    wt <- coocc_wt2(x.conts, x.cat, x.disc, bins = bins, bin_method = bin_method)
    newdf <- sweep(x.conts, 2, wt^2, function(x, y) x * sqrt(y))
    d.cont <- as.matrix(cluster::daisy(newdf, metric = "euclidean"))^2

  } else if (type == "Norm") {
    x.conts <- scale(x.cont)
    wt <- coocc_wt2(x.conts, x.cat, x.disc, bins = bins, bin_method = bin_method)
    newdf <- sweep(x.conts, 2, wt^2, function(x, y) x * sqrt(y))
    d.cont <- as.matrix(cluster::daisy(newdf, metric = "euclidean"))^2

  } else if (type == "None") {
    wt <- coocc_wt2(x.cont, x.cat, x.disc, bins = bins, bin_method = bin_method)
    newdf <- sweep(x.cont, 2, wt^2, function(x, y) x * sqrt(y))
    d.cont <- as.matrix(cluster::daisy(newdf, metric = "euclidean"))^2}

  d.Ahmad <- d.cat + d.cont

  return(list(dist=sqrt(d.Ahmad), cont=d.cont, cat=d.cat))}
