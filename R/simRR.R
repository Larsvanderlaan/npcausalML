
# Sims CATE
## no positivity
## easy RR
sim.RR <- function(n, hard = TRUE, positivity = TRUE) {
  if(!positivity & !hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    LRR <- -1 + W
    EY0W <- plogis(-1+W + sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rpois(n, EY0W * exp(A* LRR))
  }


  ## hard RR
  if(!positivity & hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    LRR <- -1 + W + sin(5*W)
    EY0W <- plogis(-1+W + sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rpois(n, EY0W * exp(A* LRR))
  }
  ##  positivity
  ## easy RR
  if(positivity & !hard) {

    W <- runif(n, -1 , 1)
    pA1W <- plogis(3*W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    LRR <- -1 + W
    EY0W <- plogis(-1+W + sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rpois(n, EY0W * exp(A* LRR))
  }
  ## hard RR
  if(positivity & hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(3*W)
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    LRR <- -1 + W + sin(5*W)
    EY0W <- plogis(-1+W + sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rpois(n, EY0W * exp(A* LRR))
  }
  return(data.table(W, A, Y, pA1W, EY1W, EY0W))
}



