library(glmnet)
library(survival)
library(MASS)

library(iterators)
library(parallel)
library(Rmpi)
library(foreach)
library(doMPI)
library(doParallel)
library(nleqslv)

cl <- startMPIcluster()
registerDoMPI(cl)

coxhdi_colcumsums <- function(x) {
  if (requireNamespace("matrixStats", quietly = TRUE)) {
    matrixStats::colCumsums(x)
  } else {
    matrix(apply(x, 2L, cumsum), nrow = nrow(x), ncol = ncol(x))
  }
}

coxhdi_logS0 <- function(eta, weight_ord) {
  eta_shift <- max(eta)
  S0 <- cumsum(weight_ord * exp(eta - eta_shift))
  if (any(!is.finite(S0)) || any(S0 <= 0)) {
    stop("Non-finite or zero risk sum; inspect extreme linear predictors.")
  }
  log(S0) + eta_shift
}

coxhdi_prepare <- function(x_cc, time_cc, status_cc, weight_cc, d, n0, n, p0,
                           obs_ind, delta_cc, beta_hat, foldid = NULL) {
  n_cc <- length(time_cc)
  beta_hat <- as.numeric(beta_hat)
  stopifnot(is.matrix(x_cc), is.numeric(x_cc), n_cc >= 2L, d >= 3L,
            nrow(x_cc) == n_cc, ncol(x_cc) == d, length(beta_hat) == d,
            length(status_cc) == n_cc, length(weight_cc) == n_cc,
            length(delta_cc) == n_cc, n0 > 0, n0 <= n,
            isTRUE(all.equal(p0, n0 / n)), all(is.finite(beta_hat)),
            all(is.finite(time_cc)), all(is.finite(weight_cc)),
            all(weight_cc > 0), all(status_cc %in% c(0, 1)),
            all(delta_cc %in% c(0, 1)))
  obs_rows <- if (is.logical(obs_ind)) which(obs_ind) else as.integer(obs_ind)
  stopifnot(length(obs_rows) == n_cc, !anyDuplicated(obs_rows),
            all(obs_rows >= 1L & obs_rows <= n))
  if (!is.null(foldid)) {
    stopifnot(length(foldid) == n, !anyNA(foldid),
              all(foldid == as.integer(foldid)), min(foldid) == 1L,
              max(foldid) >= 3L,
              identical(sort(unique(as.integer(foldid))), seq_len(max(foldid))))
  }
  
  ord <- order(time_cc, decreasing = TRUE)
  x_ord <- x_cc[ord, , drop = FALSE]
  status_ord <- status_cc[ord]
  delta_ord <- delta_cc[ord]
  weight_ord <- as.numeric(weight_cc[ord])
  eta <- as.vector(x_ord %*% beta_hat)
  if (any(!is.finite(eta))) stop("Non-finite linear predictor.")
  # A common shift cancels from every risk ratio and the projected variance.
  eta_shift <- max(eta)
  exp_xb <- exp(eta - eta_shift)
  wexp <- weight_ord * exp_xb
  S0 <- cumsum(wexp)
  if (any(!is.finite(S0)) || any(S0 <= 0)) stop("Invalid risk-set denominator.")
  event_ind <- which(status_ord == 1)
  if (!length(event_ind)) stop("There are no observed events.")
  
  # score_full[, coi] and score_full[, -coi] equal the old la_full and lb_full.
  # Keep all n rows, including zero rows: dropping them changes CV/penalties.
  score_full <- matrix(0, n, d)
  block_size <- 1024L
  for (first in seq.int(1L, d, by = block_size)) {
    cols <- seq.int(first, min(d, first + block_size - 1L))
    x_block <- x_ord[, cols, drop = FALSE]
    if (any(!is.finite(x_block))) stop("x_cc contains NA/Inf.")
    S1 <- coxhdi_colcumsums(x_block * wexp)
    score_full[obs_rows[ord[event_ind]], cols] <-
      -(x_block[event_ind, , drop = FALSE] -
          S1[event_ind, , drop = FALSE] / S0[event_ind])
  }
  list(x_ord = x_ord, status_ord = status_ord, delta_ord = delta_ord,
       weight_ord = weight_ord, eta = eta, eta_shift = eta_shift,
       exp_xb = exp_xb, wexp = wexp, S0 = S0, event_ind = event_ind,
       score_full = score_full, beta_hat = beta_hat, foldid = foldid,
       n = n, d = d, n0 = n0, p0 = p0)
}

coxhdi <- function(x_cc, time_cc, status_cc, weight_cc, d, n0, n, p0,
                   obs_ind, delta_cc, method1, beta_hat, coi = 1,
                   cache = NULL, foldid = NULL) {
  # The original positional arguments and single-method return value still work.
  # method1 = c("S", "W", "L") computes all three using one nuisance fit.
  stopifnot(length(coi) == 1L, coi == as.integer(coi), coi >= 1L, coi <= d,
            length(method1) >= 1L, all(method1 %in% c("S", "W", "L")))
  if (is.null(cache)) {
    cache <- coxhdi_prepare(x_cc, time_cc, status_cc, weight_cc, d, n0, n, p0,
                            obs_ind, delta_cc, beta_hat, foldid)
  }
  # Rebuild cache whenever x_cc, outcomes, weights, or beta_hat are changed.
  stopifnot(cache$n == n, cache$d == d, cache$n0 == n0, cache$p0 == p0,
            identical(as.numeric(beta_hat), cache$beta_hat))
  if (is.null(foldid)) foldid <- cache$foldid
  x_ord <- cache$x_ord
  status_ord <- cache$status_ord
  delta_ord <- cache$delta_ord
  weight_ord <- cache$weight_ord
  eta <- cache$eta
  exp_xb <- cache$exp_xb
  wexp <- cache$wexp
  S0 <- cache$S0
  event_ind <- cache$event_ind
  beta_hat <- cache$beta_hat
  n_cc <- nrow(x_ord)
  
  la_full <- cache$score_full[, coi]
  if (all(la_full == 0)) {
    what <- matrix(0, d - 1L, 1L)
  } else {
    lb_full <- cache$score_full[, -coi, drop = FALSE]
    cv.fit <- glmnet::cv.glmnet(x = lb_full, y = la_full,
                                standardize = FALSE, intercept = FALSE,
                                foldid = foldid, parallel = FALSE)
    fit <- cv.fit$glmnet.fit
    tmp <- match(cv.fit$lambda.min, fit$lambda)
    what <- matrix(as.numeric(fit$beta[, tmp]), nrow = d - 1L)
    rm(lb_full, cv.fit, fit)
  }
  
  # x_v = x_ord %*% v, computed using only nonzero nuisance coefficients.
  x_coi <- x_ord[, coi]
  x_v <- x_coi
  nz <- which(what[, 1L] != 0)
  if (length(nz)) {
    nuisance_cols <- seq_len(d)[-coi]
    x_v <- x_v - as.vector(x_ord[, nuisance_cols[nz], drop = FALSE] %*%
                             what[nz, , drop = FALSE])
  }
  S1_coi <- cumsum(wexp * x_coi)
  S1_v <- cumsum(wexp * x_v)
  cum_wexp_xxcoi <- cumsum(wexp * x_coi * x_v)
  
  # Exactly H_col[coi] - sum(what * H_col[-coi]), without a full H_col.
  var_H <- sum((cum_wexp_xxcoi / S0 -
                  (S1_coi / S0) * (S1_v / S0))[event_ind]) / n
  U_hat <- -sum((x_v - S1_v / S0)[event_ind]) / n
  
  var_sigma <- 0
  if (p0 < 1) {
    ord_inc <- rev(seq_len(n_cc))
    S0_inc <- S0[ord_inc]
    S1_inc <- S1_v[ord_inc]
    x_inc <- x_v[ord_inc]
    delta_inc <- delta_ord[ord_inc]
    status_inc <- status_ord[ord_inc]
    exp_xb_inc <- exp_xb[ord_inc]
    invS0 <- numeric(n_cc)
    S1_over_S0sq <- numeric(n_cc)
    is_event <- status_inc == 1
    invS0[is_event] <- 1 / S0_inc[is_event]
    S1_over_S0sq[is_event] <-
      (S1_inc[is_event] / S0_inc[is_event]) / S0_inc[is_event]
    cum_invS0 <- cumsum(invS0)
    cum_S1_over_S0sq <- cumsum(S1_over_S0sq)
    J <- which(delta_inc == 1 & status_inc == 0)
    tmp_vec <- exp_xb_inc[J] *
      (x_inc[J] * cum_invS0[J] - cum_S1_over_S0sq[J])
    # Scalars here: v' Sigma_0 v and v' Sigma_1 v.
    # The denominator is n0, NOT length(J).
    Sigma_0 <- sum(tmp_vec^2) / n0
    Sigma_1 <- (sum(tmp_vec) / n0)^2
    var_sigma <- Sigma_0 - Sigma_1
  }
  var <- var_H + (1 - p0) / p0 * var_sigma
  if (!is.finite(var) || !is.finite(var_H)) stop("Non-finite variance.")
  var <- max(as.numeric(var), 1e-6)
  var_H_safe <- max(as.numeric(var_H), 1e-6)
  atilde <- beta_hat[coi] - U_hat / var_H_safe
  Pval <- setNames(numeric(length(method1)), method1)
  
  if (any(method1 %in% c("S", "L"))) {
    eta_null <- eta - x_coi * beta_hat[coi]
  }
  if ("S" %in% method1) {
    if (beta_hat[coi] == 0) {
      S <- U_hat
    } else {
      wexp_null <- weight_ord * exp(eta_null - max(eta_null))
      S0_null <- cumsum(wexp_null)
      if (any(S0_null <= 0)) stop("Zero null risk-set denominator.")
      S1_null <- cumsum(wexp_null * x_v)
      S <- -sum((x_v - S1_null / S0_null)[event_ind]) / n
    }
    tmp <- sqrt(n) * S / sqrt(var)
    Pval[method1 == "S"] <- 2 * pnorm(-abs(tmp))
  }
  if ("W" %in% method1) {
    scale_factor <- max(var_H_safe^2 / var, 1e-8)
    tmp <- n * atilde^2 * scale_factor
    Pval[method1 == "W"] <- pchisq(tmp, df = 1, lower.tail = FALSE)
  }
  if ("L" %in% method1) {
    # Same beta1 path as the original code: eta1 = eta_null + atilde * x_v.
    eta0 <- eta_null
    eta1 <- eta0 + atilde * x_v
    logS0_1 <- coxhdi_logS0(eta1, weight_ord)
    logS0_0 <- coxhdi_logS0(eta0, weight_ord)
    tmp <- 2 * sum(((eta1 - eta0) - (logS0_1 - logS0_0))[event_ind]) *
      var_H_safe / var
    Pval[method1 == "L"] <- pchisq(tmp, df = 1, lower.tail = FALSE)
  }
  if (length(method1) == 1L) unname(Pval) else Pval
}

pval_func <- function(dim_com, x_cc = NULL, time_cc = NULL, status_cc = NULL,
                      weight_cc = NULL, obs_ind = NULL, delta_cc = NULL,
                      n, d, n0, p0, coxhdi, method1, beta_hat, cache = NULL) {
  coxhdi(x_cc, time_cc, status_cc, weight_cc, d, n0, n, p0,
         obs_ind, delta_cc, method1, beta_hat, coi = dim_com, cache = cache)
}

load("XMETABRIC.RData")
load("YMETABRIC.RData")
n <- nrow(X_dlbcl)
d <- ncol(X_dlbcl)
time_full <- Y_dlbcl[,2]
status_full <- Y_dlbcl[,1]
x_full <- X_dlbcl

n0 <- n 
p0 <- n0/n
subcohort_ind <- rep(FALSE, n)
subcohort_ind[sample(1:n, n0, replace = FALSE)] <- TRUE
saveRDS(subcohort_ind, file = sprintf("n0_%d.rds", n0))
case_ind <- (status_full == 1) 
obs_ind <- subcohort_ind | case_ind 
x_cc <- x_full[obs_ind, ]
time_cc <- time_full[obs_ind] 
status_cc <- status_full[obs_ind] 
delta_cc <- subcohort_ind[obs_ind] 
weight_full <- status_full + (1 - status_full) * (subcohort_ind / p0) 
weight_cc <- weight_full[obs_ind] 
lambdas = seq(from = 0.1, to = 3.0, by = 0.03)
lambdas = 2.05 * sqrt(log(d)) / n * sqrt(log(n) / log(n0)) * lambdas
cv.fit <- cv.glmnet(x_cc, Surv(time_cc, status_cc), family  = "cox", weights = weight_cc,          # <<< 关键修改
                    lambda  = lambdas, maxit   = 10000, nfolds  = 5, alpha   = 1)
fit     = cv.fit$glmnet.fit
tmp     = which(fit$lambda == cv.fit$lambda.min)
beta_hat    = fit$beta[,tmp]
file_name <- sprintf("beta_hat_%d.csv",n0)
write.csv(beta_hat, file = file_name, row.names = FALSE)

foldid <- sample(rep(seq_len(10L), length.out = n))
cache <- coxhdi_prepare(x_cc, time_cc, status_cc, weight_cc, d, n0, n, p0,
                        obs_ind, delta_cc, beta_hat, foldid = foldid)
method1 <- "S"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)
method1 <- "W"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)
method1 <- "L"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)


n0 <- 495
p0 <- n0/n
subcohort_ind <- rep(FALSE, n)
subcohort_ind[sample(1:n, n0, replace = FALSE)] <- TRUE
saveRDS(subcohort_ind, file = sprintf("n0_%d.rds", n0))
# subcohort_ind <- readRDS("n0_495.rds")
case_ind <- (status_full == 1) 
obs_ind <- subcohort_ind | case_ind 
x_cc <- x_full[obs_ind, ]
time_cc <- time_full[obs_ind] 
status_cc <- status_full[obs_ind] 
delta_cc <- subcohort_ind[obs_ind] 
weight_full <- status_full + (1 - status_full) * (subcohort_ind / p0) 
weight_cc <- weight_full[obs_ind] 
lambdas = seq(from = 0.1, to = 3.0, by = 0.03)
lambdas = 2.05 * sqrt(log(d)) / n * sqrt(log(n) / log(n0)) * lambdas
cv.fit <- cv.glmnet(x_cc, Surv(time_cc, status_cc), family  = "cox", weights = weight_cc,          # <<< 关键修改
                    lambda  = lambdas, maxit   = 10000, nfolds  = 5, alpha   = 1)
fit     = cv.fit$glmnet.fit
tmp     = which(fit$lambda == cv.fit$lambda.min)
beta_hat    = fit$beta[,tmp]
file_name <- sprintf("beta_hat_%d.csv",n0)
write.csv(beta_hat, file = file_name, row.names = FALSE)

foldid <- sample(rep(seq_len(10L), length.out = n))
cache <- coxhdi_prepare(x_cc, time_cc, status_cc, weight_cc, d, n0, n, p0,
                        obs_ind, delta_cc, beta_hat, foldid = foldid)
method1 <- "S"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)
method1 <- "W"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)
method1 <- "L"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)

n0 <- 990
p0 <- n0/n
subcohort_ind <- rep(FALSE, n)
subcohort_ind[sample(1:n, n0, replace = FALSE)] <- TRUE
saveRDS(subcohort_ind, file = sprintf("n0_%d.rds", n0))
# subcohort_ind <- readRDS("n0_990.rds")
case_ind <- (status_full == 1) 
obs_ind <- subcohort_ind | case_ind 
x_cc <- x_full[obs_ind, ]
time_cc <- time_full[obs_ind] 
status_cc <- status_full[obs_ind] 
delta_cc <- subcohort_ind[obs_ind] 
weight_full <- status_full + (1 - status_full) * (subcohort_ind / p0) 
weight_cc <- weight_full[obs_ind] 
lambdas = seq(from = 0.1, to = 3.0, by = 0.03)
lambdas = 2.05 * sqrt(log(d)) / n * sqrt(log(n) / log(n0)) * lambdas
cv.fit <- cv.glmnet(x_cc, Surv(time_cc, status_cc), family  = "cox", weights = weight_cc,          # <<< 关键修改
                    lambda  = lambdas, maxit   = 10000, nfolds  = 5, alpha   = 1)
fit     = cv.fit$glmnet.fit
tmp     = which(fit$lambda == cv.fit$lambda.min)
beta_hat    = fit$beta[,tmp]
file_name <- sprintf("beta_hat_%d.csv",n0)
write.csv(beta_hat, file = file_name, row.names = FALSE)

foldid <- sample(rep(seq_len(10L), length.out = n))
cache <- coxhdi_prepare(x_cc, time_cc, status_cc, weight_cc, d, n0, n, p0,
                        obs_ind, delta_cc, beta_hat, foldid = foldid)
method1 <- "S"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)
method1 <- "W"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)
method1 <- "L"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)

n0 <- 1484
p0 <- n0/n
subcohort_ind <- rep(FALSE, n)
subcohort_ind[sample(1:n, n0, replace = FALSE)] <- TRUE
saveRDS(subcohort_ind, file = sprintf("n0_%d.rds", n0))
# subcohort_ind <- readRDS("n0_1484.rds")
case_ind <- (status_full == 1) 
obs_ind <- subcohort_ind | case_ind 
x_cc <- x_full[obs_ind, ]
time_cc <- time_full[obs_ind] 
status_cc <- status_full[obs_ind] 
delta_cc <- subcohort_ind[obs_ind] 
weight_full <- status_full + (1 - status_full) * (subcohort_ind / p0) 
weight_cc <- weight_full[obs_ind] 
lambdas = seq(from = 0.1, to = 3.0, by = 0.03)
lambdas = 2.05 * sqrt(log(d)) / n * sqrt(log(n) / log(n0)) * lambdas
cv.fit <- cv.glmnet(x_cc, Surv(time_cc, status_cc), family  = "cox", weights = weight_cc,          # <<< 关键修改
                    lambda  = lambdas, maxit   = 10000, nfolds  = 5, alpha   = 1)
fit     = cv.fit$glmnet.fit
tmp     = which(fit$lambda == cv.fit$lambda.min)
beta_hat    = fit$beta[,tmp]
file_name <- sprintf("beta_hat_%d.csv",n0)
write.csv(beta_hat, file = file_name, row.names = FALSE)

foldid <- sample(rep(seq_len(10L), length.out = n))
cache <- coxhdi_prepare(x_cc, time_cc, status_cc, weight_cc, d, n0, n, p0,
                        obs_ind, delta_cc, beta_hat, foldid = foldid)
method1 <- "S"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)
method1 <- "W"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)
method1 <- "L"
pkgs <- c("MASS","glmnet","survival")
result <- foreach(coi = 1:d, .combine = cbind, .multicombine = TRUE, .packages = pkgs) %dopar% {
  pval_func(coi, x_cc,time_cc,status_cc,weight_cc, obs_ind,delta_cc, n, d, n0, p0, coxhdi, method1,beta_hat,cache)
}
file_name <- sprintf("n0%d_%s.csv",n0,method1)
write.csv(result, file = file_name, row.names = FALSE)

closeCluster(cl)
mpi.quit()