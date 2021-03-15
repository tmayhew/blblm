#' @import purrr
#' @import stats
#' @import furrr
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
#'
#'
#'
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))

#' Bag of Little Bootstraps for Linear Regression
#'
#' blblm takes in arguments for traditional linear regression and three additional arguments,
#' then returns a fitted least squares model and confidence intervals corresponding to bootstrap method.
#'
#'
#' @param formula an object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which lm is called.
#' @param m an integer, the number of equal-size groups for bag of little bootstraps algorithm.
#' @param B an integer, the number of iterations over bag of little bootstraps algorithm.
#' @param parallel logical, indicates parallelization (TRUE) or not (FALSE).
#' @param workers an integer, the number of cores to use in bag of little bootstraps algorithm.
#'
#' @export
blblm <- function(formula, data, m = 10, B = 5000, parallel = F, workers = 1){
  data_list <- split_data(data, m)
  if (parallel){
    plan(multiprocess, workers = workers)
    estimates <- future_map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B),
      .options = furrr_options(seed = TRUE))
    res <- list(estimates = estimates, formula = formula)
    class(res) <- "blblm"
    invisible(res)
  } else{
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
    res <- list(estimates = estimates, formula = formula)
    class(res) <- "blblm"
    invisible(res)
  }
}

#' Split Data
#'
#' Splits data into m parts of approximated equal sizes.
#'
#' @param data data frame with predictors and response
#' @param m integer, number of equal-size groups for bag of little bootstraps algorithm
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
lm_each_subsample <- function(formula, data, n, B) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wrong variable from the global scope.
  environment(formula) <- environment()
  m <- model.frame(formula, data)
  X <- model.matrix(formula, m)
  y <- model.response(m)
  replicate(B, lm1(X, y, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
lm1 <- function(X, y, n) {
  freqs <- as.vector(rmultinom(1, n, rep(1, nrow(X))))
  fit <- lm.wfit(X, y, freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
blbcoef <- function(fit) {
  coef(fit)
}

#' compute sigma from fit
blbsigma <- function(fit) {
  p <- fit$rank
  e <- fit$residuals
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#' Printing blblm
#'
#' Prints the formula from the blblm function.
#'
#' @param x model fit produced by blblm(formula, data, ...)
#'
#' @param ...
#'
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}

#' Sigma Value for blblm
#'
#' Outputs the value of sigma for the fitted object (after blblm algorithm)
#'
#' @param object model fit produced by blblm(formula, data, ...)
#'
#' @param confidence logical, TRUE gives a confidence interval, FALSE gives a single estimate.
#' @param level numeric, confidence interval level (only if confidence = TRUE).
#' @param ...
#'
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - level
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' Coefficients for blblm
#'
#' Outputs the coefficients for the fitted object (after blblm algorithm)
#'
#' @param object model fit produced by blblm(formula, data, ...)
#'
#' @param ...
#'
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#' Confidence Interval
#'
#' Constructs a confidence interval for parameters of the model.
#'
#' @param object model fit produced by blblm(formula, data, ...)
#'
#' @param parm character vector, parameters to give confidence interval for. Defaults to all predictors.
#' @param level numeric, confidence interval level.
#' @param ...
#'
#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(object$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' blblm model prediction
#'
#' Makes a prediction using the fitted model in object on new data with similar dimensions
#'
#' @param object model fit produced by blblm(formula, data, ...)
#'
#' @param new_data data frame with same column dimensions as in object.
#' @param confidence logical, TRUE gives a confidence interval, FALSE gives a single estimate.
#' @param level numeric, confidence interval level (only if confidence = TRUE).
#' @param ...
#'
#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}


mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}
