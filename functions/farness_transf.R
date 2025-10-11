#######################################
# Farness
#######################################
Farness <- function(dist) {

  indnz <- which(dist > 1e-10)
  farnz <- dist[indnz]
  farloc <- median(farnz, na.rm = TRUE)
  farsca <- mad(farnz, na.rm = TRUE)

  if (farsca < 1e-10) {farsca <- sd(farnz, na.rm = TRUE)}
  sfar <- scale(farnz, center = farloc, scale = farsca)
  YJout  <- cellWise::transfo(X = sfar, type = "YJ", robust = TRUE, standardize = FALSE, checkPars = list(silent = TRUE))
  xt <- YJout$Y

  tfarloc <- median(xt, na.rm = TRUE)
  tfarsca <- mad(xt, na.rm = TRUE)
  zt <- scale(xt, tfarloc, tfarsca)
  probs <- rep(0, length(dist))
  probs[indnz] <- pnorm(zt)

  return(farness = probs)}
