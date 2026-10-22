library(glmnet)
library(survival)
library(mvtnorm)
library(glasso)
library(lpSolve)
library(MASS)
# rm(list = ls()); gc()

library(iterators)
library(parallel)
library(Rmpi)
library(foreach)
library(doMPI)
library(doParallel)
library(nleqslv)
cl <- startMPIcluster()
registerDoMPI(cl)

coxhdi <- function(x_cc,time_cc,status_cc,weight_cc,d,n0,n,p0,obs_ind,delta_cc,method1,coi=1)
{
  Pval = c()
  lambdas = seq(from = 0.1, to = 1.5, by = 0.03)
  lambdas = 2.05 * sqrt(log(d)) / n * sqrt(log(n) / log(n0)) * lambdas
  cv.fit <- cv.glmnet(x_cc, Surv(time_cc, status_cc), family  = "cox", weights = weight_cc,          # <<< 关键修改
                      lambda  = lambdas, maxit   = 10000, nfolds  = 5, alpha   = 1)
  fit     = cv.fit$glmnet.fit
  tmp     = which(fit$lambda == cv.fit$lambda.min)
  beta_hat    = fit$beta[,tmp]
  
  n_cc = length(time_cc)
  Vs  = matrix(0, d, d)
  Hs  = Vs
  la = rep(0, n_cc)
  lb = matrix(0, n_cc, d-1)
  i = 1
  while (i <= n_cc){
    if (status_cc[i]==1){
      ind = which(time_cc >= time_cc[i])
      S0 = 0
      S1 = rep(0, d)
      S2 = matrix(0, d, d)
      for (j in ind){
        wj  = weight_cc[j]
        exb = exp(x_cc[j, ] %*% beta_hat)
        tmp = wj * exb
        S0 = S0 + tmp
        S1 = S1 + tmp %*% t(x_cc[j, ])
        S2 = S2 + as.numeric(tmp) * (x_cc[j, ] %*% t(x_cc[j, ]))
      }
      la[i] = -(x_cc[i, coi] - S1[coi] / S0) 
      lb[i, ] = -(x_cc[i, -coi] - S1[-coi] / as.numeric(S0))
      V  = as.numeric(S0) * S2 - t(S1) %*% S1
      Hs = Hs + V / (n * as.numeric(S0)^2)
    }
    i = i + 1
  }
  
  la_full <- rep(0, n)
  lb_full <- matrix(0, n, d-1)
  la_full[obs_ind] = la
  lb_full[obs_ind,] = lb
  
  cv.fit = cv.glmnet(x = lb_full, y = la_full, standardize = FALSE, intercept  = FALSE)
  fit = cv.fit$glmnet.fit
  tmp = which(fit$lambda == cv.fit$lambda.min)
  what = fit$beta[, tmp]
  what = matrix(what,nrow=d-1)
  
  Sigma_0 = matrix(0, d, d)
  Sigma_1_vec = rep(0, d)
  J = which(delta_cc==1 & status_cc==0)
  for (j in J){
    sum_i_vec = rep(0, d)
    ind = which(time_cc <= time_cc[j])
    for (i in ind){
      if (status_cc[i]==1){
        ind_l = which(time_cc >= time_cc[i])
        S0 = 0
        S1 = rep(0, d)
        for (l in ind_l){
          wl = weight_cc[l]
          exb = exp(x_cc[l, ] %*% beta_hat)
          tmp = wl * exb
          S0 = S0 + tmp
          S1 = S1 + tmp %*% t(x_cc[l, ])
        }
        S0 = as.numeric(S0)
        sum_i_vec = sum_i_vec + (x_cc[j, ] - S1/S0) / S0
      }
    }
    if (all(sum_i_vec == 0)) next
    exb = exp(x_cc[j, ] %*% beta_hat)
    tmp_vec = matrix(as.numeric(exb) * sum_i_vec, ncol = 1)
    Sigma_0 = Sigma_0 + (tmp_vec %*% t(tmp_vec)) / n0
    Sigma_1_vec = Sigma_1_vec + tmp_vec / n0
  }
  Sigma_1 = Sigma_1_vec %*% t(Sigma_1_vec)
  v <- c(1, -as.numeric(what))
  var_sigma = t(v) %*% (Sigma_0-Sigma_1) %*% v
  
  var_H = Hs[coi,coi] - t(what)%*%Hs[-coi, coi]
  var = var_H + (1-p0)/p0 * var_sigma
  var = max(as.numeric(var), 1e-6)
  
  
  if (method1 == "S"){
    # Decorrelated Score
    beta_null = beta_hat
    beta_null[coi] = 0
    la_null = rep(0, n_cc)
    lb_null = matrix(0, n_cc, d-1)
    i = 1
    while (i <= n_cc){
      if (status_cc[i]==1){
        ind = which(time_cc >= time_cc[i])
        S0 = 0
        S1 = rep(0, d)
        S2 = matrix(0, d, d)
        for (j in ind){
          wj  = weight_cc[j]
          exb = exp(x_cc[j, ] %*% beta_null)
          tmp = wj * exb
          S0 = S0 + tmp
          S1 = S1 + tmp %*% t(x_cc[j, ])
        }
        la_null[i] = -(x_cc[i, coi] - S1[coi] / S0) 
        lb_null[i, ] = -(x_cc[i, -coi] - S1[-coi] / as.numeric(S0))
      }
      i = i + 1
    }
    lbt = colSums(lb_null)
    lbt = matrix(lbt,nrow = d-1)
    S = sum(la_null) - t(what)%*%(lbt)
    S = S/n
    tmp  = sqrt(n)*S/sqrt(var)
    pval__  = 2*pnorm(-abs(tmp))
  } else if (method1 == "W") {
    # Decorrelated Wald
    lbt = colSums(lb)
    lbt = matrix(lbt,nrow = d-1)
    U_hat = sum(la) - t(what)%*%(lbt)
    U_hat = U_hat/n
    var_H = max(as.numeric(var_H), 1e-6)
    atilde = beta_hat[coi] - U_hat/var_H
    tmp = (n)*atilde^2*(max(var_H^2/var, 1e-8))
    pval__  = 1-pchisq(tmp,1)
  } else if (method1 == "L") {
    # Deccorelated LHR
    lbt = colSums(lb)
    lbt = matrix(lbt,nrow = d-1)
    U_hat = sum(la) - t(what)%*%(lbt)
    U_hat = U_hat/n
    var_H = max(as.numeric(var_H), 1e-6)
    atilde = beta_hat[coi] - U_hat/var_H
    betal = c(atilde, beta_hat[-coi]-as.numeric(atilde)*what)
    if (coi == 1)
    {
      beta1 = c(atilde,beta_hat[-coi]-as.numeric(atilde)*what)
    } else if (coi == d) {
      beta1 = c(beta_hat[-coi]-as.numeric(atilde)*what,atilde)
    } else {
      beta1 = c(beta_hat[1:(coi-1)]-as.numeric(atilde)*what[1:(coi-1)],atilde,beta_hat[(coi+1):d]-as.numeric(atilde)*what[(coi):(d-1)])
    }
    
    lossa  = 0    #Log-Likelihood under alternative
    i = 1
    while( i<=n_cc)
    {
      if (status_cc[i]==1){
        ind = which(time_cc >= time_cc[i])
        S0 = 0
        for (j in ind){
          tmp = exp(x_cc[j,]%*%beta1)
          S0 = S0 + as.numeric(tmp) * weight_cc[j]
        }
        lossa = lossa + 1/n * (x_cc[i,]%*%beta1 - log(S0))
      }
      i = i + 1
    }
    
    lossn  = 0  #Log-Likelihood under null
    beta_null = beta_hat
    beta_null[coi] = 0
    i = 1
    while(i<=n_cc){
      if (status_cc[i]==1){
        ind = which(time_cc >= time_cc[i])
        S0 = 0
        for (j in ind){
          tmp = exp(x_cc[j,]%*%beta_null)
          S0 = S0 + as.numeric(tmp) * weight_cc[j]
        }
        lossn = lossn + 1/n * (x_cc[i,]%*%beta_null - log(S0))
      }
      i = i + 1
    }
    tmp   = n*lossa - n*lossn
    tmp   = 2*tmp*var_H/var
    pval__  = 1-pchisq(tmp,1)
  }
  Pval = c(pval__)
}

pval_func <- function(iter, n, d, s, Sigma, n0, p0, coxhdi, method1, rej_data){ 
  beta_star <- rep(0, d)
  active_set <- sample(2:d, s)
  beta_star[active_set] <- 1
  beta_star[1] <- rej_data
  x_full <- MASS::mvrnorm(n = n, mu = rep(0, d), Sigma = Sigma) 
  linpred <- as.vector(x_full %*% beta_star) 
  T_true <- rexp(n, rate = exp(linpred)) 
  U <- runif(n, 1, 3) 
  C <- rexp(n, rate = 1 / (U * exp(linpred))) 
  time_full <- pmin(T_true, C) 
  status_full <- as.numeric(T_true <= C) 
  subcohort_ind <- rep(FALSE, n) 
  subcohort_ind[sample(1:n, n0, replace = FALSE)] <- TRUE 
  case_ind <- (status_full == 1) 
  obs_ind <- subcohort_ind | case_ind 
  x_cc <- x_full[obs_ind, ] 
  time_cc <- time_full[obs_ind] 
  status_cc <- status_full[obs_ind] 
  delta_cc <- subcohort_ind[obs_ind] 
  weight_full <- status_full + (1 - status_full) * (subcohort_ind / p0) 
  weight_cc <- weight_full[obs_ind] 
  pval_ <- coxhdi(x_cc,time_cc,status_cc,weight_cc,d,n0,n,p0,obs_ind,delta_cc,method1,coi=1) 
  return(pval_) } 


d <- 500
n <- 200
s <- 2
n0 <- 50
p0 <- n0 / n
rho <- 0.4
iters <- 1000
Sigma <- outer(1:d, 1:d, function(j,k) rho^abs(j-k))
method1_list = c("S", "W", "L")
rej_data_list = seq(0, 0.55, by = 0.01)
for (rej_data in rej_data_list){
  for (method1 in method1_list){
    result <- foreach(iter = 1:iters, .combine = cbind, .multicombine = TRUE,
                    .packages = c("MASS","glmnet","survival","mvtnorm","glasso","lpSolve")) %dopar% {
                    pval_func(iter, n, d, s, Sigma, n0, p0, coxhdi, method1, rej_data)
                    }
    file_name <- sprintf("rej_d%d_n0%03d_rho%02d_Dirac_%s_s%d_%02d.csv",
                          d,n0,as.integer(round(rho * 100)),method1,s,round(rej_data*100))
    write.csv(result, file = file_name, row.names = FALSE)
  }
}

closeCluster(cl) 
mpi.quit()



#squeue -u $USER
#scancel -u $USER
#sbatch test300.sh
##SBATCH --mem=8G
#SBATCH --mem-per-cpu=4G        # 每核分配内存
