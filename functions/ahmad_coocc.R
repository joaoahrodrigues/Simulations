#######################################
# coocc_dis
#######################################
coocc_dis<-function(x.cat, x.cont, x.disc, only.cat = FALSE, js = FALSE, bins = NULL, bin_method = NULL) {

  bin_method <- match.arg(bin_method, choices = c("equalwidth", "weighted"))
  if (only.cat == TRUE) {
    df <- x.cat
  } else {
    df <- cbind(x.cat, x.disc)}

  delta <- vector("list", ncol(df))
  for (i in 1:ncol(df)){
    delta[[i]] <- 0
    for (j in 1:(ncol(df)-1)){
      tmp <- xtabs(~ df[,i] + df[,-i,drop=FALSE][,j])
      tmp <- tmp/rowSums(tmp)
      if (js) { # With Jensen Shannon
        js_mat <- outer(1:nrow(tmp), 1:nrow(tmp), Vectorize(function(r1, r2)
          as.matrix(distance(rbind(tmp[r1, ], tmp[r2, ]), method = "jensen-shannon"))))

        delta[[i]] <- delta[[i]] + js_mat}

      else { # Without Jensen Shannon
        delta[[i]] <- delta[[i]] + as.matrix(dist(tmp, "manhattan"))}}
    delta[[i]] <- 0.5*delta[[i]]/(ncol(df) - 1)

    if (js) {
      dimnames(delta[[i]]) <- list(levels(as.factor(df[,i])), levels(as.factor(df[,i])))}}

  dist_mat <- matrix(0, nrow(df), nrow(df))

  for (k in 1:ncol(df)){
    tmp_col <- df[,k]

    for (i in 1:(nrow(df) - 1)) {
      for (j in (i + 1):nrow(df)) {
        val <- delta[[k]][as.character(tmp_col[i]), as.character(tmp_col[j])] ^ 2
        dist_mat[i, j] <- dist_mat[i, j] + val
        dist_mat[j, i] <- dist_mat[i, j]}}}

  diag(dist_mat) <- 0

  return(dist_mat)}
