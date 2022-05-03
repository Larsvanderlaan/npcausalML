# CATE
#' @export
outcome_function_plugin_CATE <- function(A, Y, EY1W, EY0W, pA1W) {
  EY1W - EY0W
}
#' @export
weight_function_plugin_CATE <- function(A, Y, EY1W, EY0W, pA1W) {
  rep(1, length(A))
}
#' @export
outcome_function_IPW_CATE <- function(A, Y, EY1W, EY0W, pA1W) {
  pA0W <- 1- pA1W
  Y * (A/pA1W - (1-A)/(pA0W))
}
#' @export
weight_function_IPW_CATE <- function(A, Y, EY1W, EY0W, pA1W) {
  #pA <- ifelse(A==1, pA1W, 1 - pA1W)
  #1 / pA
  return(rep(1,length(A)))
}


#' @export
design_function_sieve_plugin_CATE <- function(X,A , Y, EY1W , EY0W , pA1W ) {
   (A*X-(1-A)*X)/ifelse(A==1,pA1W, 1- pA1W)
}
#' @export
design_function_sieve_IPW_CATE <- function(X,A , Y, EY1W , EY0W , pA1W ){
  pA0W <- 1-pA1W
  cbind(EY1W/pA1W * X, X* EY0W/pA0W)
}
#' @export
weight_function_sieve_plugin_CATE <- function(A , Y, EY1W , EY0W , pA1W ){
  rep(1, length(A))#1/ifelse(A==1,pA1W, 1- pA1W)
}
#' @export
weight_function_sieve_IPW_CATE <- function(A , Y, EY1W , EY0W , pA1W ){
  return(rep(1,length(A)))
}


#' @export
efficient_loss_function_CATE <- function(V, theta, A , Y, EY1W , EY0W , pA1W, oracle = FALSE){
  pA <-  ifelse(A==1,pA1W, 1- pA1W)
  EY <- ifelse(A==1, EY1W, EY0W)
  CATE <- EY1W - EY0W
  loss <- theta^2 - 2*theta*(CATE + (1/pA)*(A - (1-A))*(Y - EY))


  return(loss)
}


## RR
#' @export
outcome_function_plugin_LRR <- function(A, Y, EY1W, EY0W, pA1W) {
  #print("outcome plugin")


  #print(quantile((EY1W  / (EY1W + EY0W))))
  EY1W  / (EY1W + EY0W)
}
#' @export
weight_function_plugin_LRR <- function(A, Y, EY1W, EY0W, pA1W) {
  (EY1W + EY0W)
  #print("weight")
  #print(quantile((  (EY1W + EY0W))))
  (EY1W + EY0W)
}
#' @export
outcome_function_IPW_LRR <- function(A, Y, EY1W, EY0W, pA1W) {
  A
}
#' @export
weight_function_IPW_LRR <- function(A, Y, EY1W, EY0W, pA1W) {
  pA <- ifelse(A==1, pA1W, 1 - pA1W)
  Y / pA
}

#' @export
design_function_sieve_plugin_LRR <- function(X,A , Y, EY1W , EY0W , pA1W ) {
  cbind(A* X, X*(1-A))
}
#' @export
design_function_sieve_IPW_LRR <- function(X,A , Y, EY1W , EY0W , pA1W ){
  pA0 <- 1-pA1W

  cbind(EY1W/pA1W * X, EY0W/pA0 * X)
}
#' @export
weight_function_sieve_plugin_LRR <- function(A , Y, EY1W , EY0W , pA1W ){
  1/ifelse(A==1,pA1W, 1- pA1W)
}
#' @export
weight_function_sieve_IPW_LRR <- function(A , Y, EY1W , EY0W , pA1W ){
  return(rep(1,length(A)))
}

#' @export
efficient_loss_function_LRR <- function(theta, A , Y, EY1W , EY0W , pA1W){
  LRR <- theta
  EY <- ifelse(A==1, EY1W, EY0W)
  plugin_risk <- (EY0W + EY1W) * log(1 + exp(LRR)) - EY1W * LRR
  score_comp <- (A/pA1W)*(log(1 + exp(LRR)) - LRR)*(Y - EY) + ((1-A)/(1-pA1W))*(log(1 + exp(LRR)) - LRR)*(Y - EY)
  plugin_risk + score_comp
}
