
make_super_learner <- function(candidates, training_task, likelihood, loss_generator = make_eff_loss) {
  loss_function <- loss_generator(training_task, likelihood)
  Lrnr_sl$new(candidates,  metalearner = Lrnr_cv_selector$new(loss_function = loss_function))
}


