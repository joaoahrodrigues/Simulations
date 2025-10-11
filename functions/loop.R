#######################################
# LoOP function
#######################################
LOOP <- function(dist, K=5, lambda=3){
  
  if(!inherits(dist,"dist")){
    stop('dist input is not class `dist`')}
  
  n <- attr(dist,"Size")
  
  if(!is.numeric(K)){
    stop('k input must be numeric')}
  if(K>=n||K<1){
    stop('k input must be less than number of observations and greater than 0')}
  if(!is.numeric(lambda)){
    stop('lambda input must be numeric')}
  
  dist.obj <- dbscan::kNN(dist, K)
  nnSD <- apply(dist.obj$dist, 1, function(x){sqrt((sum(x^2)/K))})
  
  plof <- NULL
  
  for(i in 1:n){
    plof[i] <- (nnSD[i]/mean(nnSD[dist.obj$id[i,]]))-1}
  nplof <- lambda*sqrt(sum(plof^2)/n)
  loop <- pracma::erf(plof/(nplof*sqrt(2)))
  loop[loop<0] <- 0
  
  return(loop)}
