
# Sims CATE
## no positivity
## easy RR
sim.RR <- function(n, hard = TRUE, positivity = TRUE) {
  if(!positivity & !hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(W + sin(5*W))
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    LRR <- -1 + W
    EY0W <- plogis(-1+W + sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rbinom(n, size = 1, EY0W * (1-A) + A* EY1W)
  }


  ## hard RR
  if(!positivity & hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(W+ sin(5*W))
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    LRR <- 0.1*(2*abs(W) + -1 +W - 2*abs(W)^3 + 1.5*sin(5*W) - cos( 1 + 7*W)- 0.5*sin(9*W+1))
    EY0W <- plogis(-1+W - sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rbinom(n, size = 1, EY0W * (1-A) + A* EY1W)
  }
  ##  positivity
  ## easy RR
  if(positivity & !hard) {

    W <- runif(n, -1 , 1)
    pA1W <- plogis(4*W + sin(5*W))
    A <- rbinom(n, 1 ,  pA1W)
    quantile(plogis(W))
    LRR <- -1 + W
    EY0W <- plogis(-1+W - sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rbinom(n, size = 1, EY0W * (1-A) + A* EY1W)
  }
  ## hard RR
  if(positivity & hard) {
    W <- runif(n, -1 , 1)
    pA1W <- plogis(4*W + sin(5*W))
    A <- rbinom(n, 1 ,  pA1W)
    LRR <- 0.1*(2*abs(W) + -1 +W - 2*abs(W)^3 + 1.5*sin(5*W) - cos( 1 + 7*W)- 0.5*sin(9*W+1))
    EY0W <- plogis(-1+W - sin(5*W) + 1/(W + 1.2))
    EY1W <- EY0W * exp(LRR)
    Y <- rbinom(n, size = 1, EY0W * (1-A) + A* EY1W)
  }
  return(data.table(W, A, Y, pA1W, EY1W, EY0W))
}




