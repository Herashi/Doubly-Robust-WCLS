setwd("~/Documents/GitHub/DML-WCLS/Doubly Robust WCLS clipped bandit")
source("init.R")


## set number of Monte Carlo replicates
M <- 1000


## set number of threads to use for parallel processing and the random seed
## (nb: these two values ensure that the results are replicable)
cores <- 4
seed <- 123

cl <- makeCluster(getOption("cl.cores", cores))
clusterEvalQ(cl, source("init.R"))
registerDoParallel(cl)

## create time-varying contextual information
high_d = 20
Control_var = paste0("Control_var", 1:high_d,collapse = " + ",sep = "")

sim.omit <- function() {
  out <- NULL
  ## low, medium and high degrees of moderation by state
  for (b in 0.2) {
    ## number of independent individuals
    for (n in 3) {
      ## number of time points observed for each individual
      for (tmax in 150) {
        clusterSetRNGStream(cl, seed)
        out <-sim_wc(n, tmax, M, high_d = 20,
                             ## regress response on state and proximal treatment,
                             ## ignoring the underlying interaction between the two
                             y.formula = list(w = as.formula(paste0("y ~ state + I(a - pn) + ", Control_var))),
                             contrast_vec = c(0,0,1, rep(0,high_d)),
                             y.names = c(w = "Doubly Robust Large T"),
                             ## term labels for proximal treatment
                             y.label = list(w = "I(a - pn)"),
                             ## specify weights and working correlation structure
                             y.args = list(w = list(wn = "pn", wd = "prob")),
                             ## specify weight numerator model
                             a.formula = list(pn = a ~ 1),
                             a.names = c(pn = "Intercept-only"),
                             ## use default generative model, but with the specified
                             ## level of moderation by the time-varying state
                             beta0 = c(-0.2, 0, 0, b, 0))
      }
    }
  }
  out
}




omit <- sim.omit()
save(omit,file = "test.RData")

stopCluster(cl)