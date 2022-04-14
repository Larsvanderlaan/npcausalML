n <- 250
W <- runif(n, min = -1,  max = 1)
A <- rbinom(n, size = 1, prob = plogis(W))
Y <- rbinom(n,1,plogis(A * (1 + W + 2*W^2) + sin(5*W)))
propW <- 1/rnorm(n,0.5,0.1)
censoringW <-  runif(n, min = 1,  max = 5)
data <- data.frame(W,A,Y, propW, censoringW)
library(npRR)
# To install the sl3 library run the following line of code:
## devtools::install_github("tlverse/sl3", ref = "devel")
if(!require(sl3)) {
  devtools::install_github("tlverse/sl3", ref = "devel")
}
### fit a working log-linear model
# formula_LRR is a formula object specifying the linear form of the log relative risk.
# For the vaccine applications, this would specify a linear form for the log modified vaccine efficacy
# The following formula uses the working model log(E[Y|A=1,W]/E[Y|A=0,W]) = beta0 + beta1 W.
formula_LRR <- ~ 1 + W
# The inverse probability weights are probably 1 / {P(sampled)*P(not_censored)}
IPW_weights <- propW * censoringW

# The machine-learning algorithms for estimating E[Y|A,W] and P(A=1 |W)
# can be specified by passing in a sl3_Learner object from the package github tlverse/sl3
lrnr_stack <-  Stack$new(list( Lrnr_glmnet$new() ,
                            Lrnr_gam$new() ,
                             Lrnr_earth$new() ,
            Lrnr_xgboost$new(max_depth = 4),
            Lrnr_xgboost$new(max_depth = 5)))
# Stratifies estimation by treatment group
lrnr_stack_stratified <- Lrnr_stratified$new(lrnr_stack, "A")

lrnr_cv_stack <- Lrnr_cv$new(lrnr_stack)
lrnr_superLearner_A <- make_learner(Pipeline, lrnr_cv_stack, Lrnr_cv_selector$new(loss_loglik_binomial))
sl3_Learner_pA1W <- lrnr_superLearner

lrnr_cv_stack <- Lrnr_cv$new(lrnr_stack_stratified)
lrnr_superLearner_Y <- make_learner(Pipeline, lrnr_cv_stack, Lrnr_cv_selector$new(loss_loglik_binomial))
sl3_Learner_EYAW <- lrnr_superLearner

# Make sure that W is a named matrix
fit <- npRRWorkingModel(formula_LRR, data.frame(W=W), A, Y, weights = IPW_weights, sl3_Learner_EYAW = sl3_Learner_EYAW, sl3_Learner_pA1W = sl3_Learner_pA1W )
coef(fit)



