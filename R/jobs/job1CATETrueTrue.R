library(future)
plan(multicore, workers = 8)
ns <- c(100, 250, 500, 750, 1000, 5000, 10000)

output <- lapply(ns, function(n) {
  out <- estCATE(n,  lrnr_stack, 100, sim.CATE, list_of_sieves_uni, positivity = FALSE, hard = FALSE)
  out$n <- n
  out
})
output <- do.call(rbind, output)
save(output, "CATEsims_posFALSEhardFALSE.RData")
