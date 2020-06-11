#' @import purrr
#' @import stats
#' @import parallel
#' @aliases NULL
#' @importFrom magrittr %>%
#' @importFrom utils capture.output
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' Generate linear regression with boostrap
#'
#' @param formula regression formula
#' @param data the data you want to regression
#' @param m number of subdata
#' @param B number of bootstrips
#'
#' @return blblm
#' @export
#' @examples
#' blblm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
blblm <- function(formula, data, m = m, B = B) {
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}
#' Generate linear regression with boostrap
#'
#' @param formula regression formula
#' @param data data we are dealing with
#' @param m number of subdata
#' @param B number of bootstrips
#' @param cl cluster/integer
#'
#' @return blblm_par
#' @export
blblm_par <- function(formula, data, m = m, B = B,cl) {
  data_list <- split_data(data, m)
  estimates <- parLapply(cl,data_list,fun=lm_each_subsample,formula = formula, n = nrow(data), B = B)
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}

#' Generate linear regression with boostrap
#'
#' @param formula regression formula
#' @param data the data you want to regression
#' @param m number of subdata
#' @param B number of bootstrips
#'
#' @return blblm
#' @export
#' @examples
#' blblm_glm(mpg ~ wt * hp, data = mtcars, m = 3, B = 100)
blblm_glm <- function(formula, data, m = m, B = B) {
  data_list <- split_data(data, m)
  estimates <- map(
    data_list,
    ~ glm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}
#' Generate linear regression with boostrap
#'
#' @param formula regression formula
#' @param data data we are dealing with
#' @param m number of subdata
#' @param B number of bootstrips
#' @param cl cluster/integer
#'
#' @return blblm_par
#' @export
blblm_glm_par <- function(formula, data, m = m, B = B,cl) {
  data_list <- split_data(data, m)
  estimates <- parLapply(cl,data_list,fun=glm_each_subsample,formula = formula, n = nrow(data), B = B)
  res <- list(estimates = estimates, formula = formula)
  class(res) <- "blblm"
  invisible(res)
}


#' split data into m parts of approximated equal sizes
#' @param m number of subdata
#' @param data data we are dealing with
#' @export
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
#' @param formula regression formula
#' @param n integer
#' @param B number of bootstrips
#' @param data data we are dealing with
#' @export
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}

#' compute the regression estimates for a blb dataset
#' @param formula regression formula
#' @param data data we are dealing with
#' @param n integer
#' @export
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
#' @param formula regression formula
#' @param data data we are dealing with
#' @param freqs frequency
#' @export
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  object <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(object), sigma = blbsigma(object))
}
#' compute the estimates
#' @param formula regression formula
#' @param n integer
#' @param B number of bootstrips
#' @param data data we are dealing with
#' @export
glm_each_subsample <- function(formula, data, n, B) {
  replicate(B, glm_each_boot(formula, data, n), simplify = FALSE)
}

#' compute the regression estimates for a blb dataset
#' @param formula regression formula
#' @param data data we are dealing with
#' @param n integer
#' @export
glm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  glm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
#' @param formula regression formula
#' @param data data we are dealing with
#' @param freqs frequency
#' @export
glm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  object <- glm(formula, data, weights = freqs,family = gaussian)
  coef = blbcoef(object) #glm does not have sigma
}

#' compute the coefficients from fit
#' @param object fit model
#' @export
blbcoef <- function(object) {
  coef(object)
}


#' compute sigma from fit
#' @param object fit model
#' @export
blbsigma <- function(object) {
  p <- object$rank
  y <- model.extract(object$model, "response")
  e <- fitted(object) - y
  w <- object$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}

#print function
#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}

#function to calculate sigma
#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}
#find the coeficient
#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}

#Find the confidence interval
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
#prediction
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

#calculation lower to upper bound
mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(object = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}
#calculation of mean
map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}
#column bind
map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}
#row bind
map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}






