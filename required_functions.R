library(reshape2)
library(ggplot2)
library(rstan)
library(dplyr)
options(mc.cores = parallel::detectCores())

aij_interval <- function(point, stan_object, row, CI= .95) {
  calc_alpha <- function(mu0,  mu1, sigma, mix, point) {
    mu0_d <- mix*dnorm(point, mean = mu0, sd = sigma)
    mu1_d <- (1-mix)*dnorm(point, mean = mu1, sd = sigma) + .000001
    
    alpha <- mu1_d/(mu0_d + mu1_d)
    return(alpha)
  }
  
  alpha_wrap <- function(store_list, points, additive = FALSE) {
    if(!additive) {
      calc_alpha(store_list[["mu0"]], store_list[["mu1"]], store_list[["sigma"]],
                 store_list[["mix"]], point = points)
    } else {
      calc_alpha(store_list[["mu0"]], store_list[["mu0"]] + store_list[["mu1"]], store_list[["sigma"]],
                 store_list[["mix"]], point = points)
    }
  }
  
  
  
  fit_frame <- as.data.frame(stan_object)
  
  mu0 <- fit_frame[[paste0("mu0[", row, "]")]]
  mu1 <- fit_frame[[paste0("mu1[", row, "]")]]
  mix <- fit_frame[[paste0("theta[",row,"]")]]
  sigma <- fit_frame[[paste0("sigma[",row,"]")]]
  CI <- 1- CI
  
  alphas <- mapply(calc_alpha, mu0 = mu0, mu1 = mu1, mix = mix, sigma = sigma, point = point)
  return_vec <- alphas[rank(alphas) > length(alphas)*CI/2  & rank(alphas) < length(alphas)*(1-CI/2)]
  return(c(min(return_vec), max(return_vec)))
}
param_interval <- function(stan_object, row, CI = .95) {
  fit_frame <- as.data.frame(stan_object)
  CI <- 1-CI
  
  get_interval <- function(vec, CI) {
    vec <- vec[rank(vec) > CI*length(vec) & rank(vec) < (1-CI)*length(vec)]
    return(c(min(vec), max(vec)))
  }
  
  mu0 <- fit_frame[[paste0("mu0[", row, "]")]]
  mu1 <- fit_frame[[paste0("mu1[", row, "]")]]
  mixes <- fit_frame[[paste0("theta[",row,"]")]]
  sigma <- fit_frame[[paste0("sigma[",row,"]")]]
  
  return_list <- lapply(list("mu0" = mu0, mu1 = mu1, theta = mixes, sigma = sigma), get_interval, CI = CI)
}
rank_method <- function(sample) {
  ifelse(rank(sample)/length(sample) <.5, 1, 0)
}
get_mclust <- function(model_object) {
  params_list <- model_object$parameters
  means <- as.numeric(params_list$mean[order(abs(params_list$mean))])
  variance <- params_list$variance$sigmasq
  mix <- params_list$pro[order(abs(params_list$mean))][1]
  bic <- model_object$bic
  
  return(list(mu0 = means[1], mu1 = means[2], sigma = sqrt(variance) , mix = mix, bic = bic))
}
#parameters for stan model when mixture not fixed at zero
calc_alpha <- function(mu0,  mu1, sigma, mix, point) {
  mu0_d <- mix*dnorm(point, mean = mu0, sd = sigma)
  mu1_d <- (1-mix)*dnorm(point, mean = mu1, sd = sigma) + .000001
  
  alpha <- mu1_d/(mu0_d + mu1_d)
  return(alpha)
}
alpha_wrap <- function(store_list, points, additive = FALSE) {
  if(!additive) {
    calc_alpha(store_list[["mu0"]], store_list[["mu1"]], store_list[["sigma"]],
               store_list[["mix"]], point = points)
  } else {
    calc_alpha(store_list[["mu0"]], store_list[["mu0"]] + store_list[["mu1"]], store_list[["sigma"]],
               store_list[["mix"]], point = points)
  }
}
generate_sample <- function(n, mu0, mu1, sd, theta0) {
  mix1 <- rnorm(n, mu0, sd)
  mix2 <- rnorm(n, mu1, sd)
  sampling <- rbinom(n,1, 1 - theta0)
  return_vec <- rep(NA, n)
  return_vec[sampling == 0] <- mix1[sampling == 0]
  return_vec[sampling == 1] <- mix2[sampling == 1]
  return(list(sampling = return_vec, labels = sampling))
}
rank_method <- function(sample) {
  ifelse(rank(sample)/length(sample) <.5, 1, 0)
}

sim <- function(mu0, mu1, theta0, sd, n=170){
  require(rstan)
  require(mclust)

  extract_fix_param <- function(model_vec, row) {
    
    mu0 <- models[paste0("mu0[", row, "]")]
    mu1 <- models[paste0("mu1[", row, "]")]
    mixes <- models[paste0("theta[",row,"]")]
    sigma <- models[paste0("sigma[",row,"]")]
    
    return(list(mu0 = mu0, mu1 = mu1, sigma = sigma, mix = mixes))
  }
  
  param_matrix <- matrix(c(mu0,mu1,theta0,sd),nrow=length(mu0))
  n_mix <- sample(1:5, nrow(param_matrix), replace=T)
  
  stor_sample <- matrix(nrow=sum(n_mix), ncol=n)
  stor_labels <- matrix(nrow=sum(n_mix), ncol=n)
  current_row = 1
  return_param_matrix <- matrix(0, nrow = sum(n_mix), ncol = 4)
  
  eval_accuracy <- function(sample_points, labels, model, additive = F, falsePos = F) {
    aij <- alpha_wrap(model, sample_points, additive = additive)
    
    label <- ifelse(aij > .5, 1, 0)
    accuracy <- mean(label == labels)
    falsePositive <- sum(label != labels & label == 1)/sum(labels == 0)
    if(falsePos) {
      falsePositive
    } else {
      accuracy
    }
  }
  
  eval_entropy <- function(sample_points, labels, model, additive = F) {
    aij <- alpha_wrap(model, sample_points, additive = additive)
    labels = ifelse(aij > .5, 1, 0)
    -sum((labels*log(aij + .00000001) + (1 - labels)*log(1 - aij + .000001)))
  }
  
  for (i in 1:nrow(param_matrix)){
    for(m in 1:n_mix[i]) {
      temp <- generate_sample(n,mu0[i],mu1[i],sd[i],theta0[i])
      stor_sample[current_row, ] <- 
        temp[[1]]
      stor_labels[current_row,] <- 
        temp[[2]]
      return_param_matrix[current_row, ] <- param_matrix[i,]
      current_row = current_row+1
    }
  }
  
  stanFeed <- list(N = nrow(stor_sample), J = ncol(stor_sample), y = stor_sample)
  fit <- stan(file = "simple_multivar.stan", data= stanFeed, iter = 500, chains = 4 , control = list(max_treedepth = 15))
  fit_frame <- as.data.frame(fit)
  models <- param_vec <- apply(X = fit_frame,MARGIN = 2, mean)
  
  mclust_store <- apply(stor_sample, MARGIN = 1, Mclust, G = 2, modelNames = "E")
  
  #calculate performance statistics
  m_m_perf <- sapply(1:nrow(stor_labels), function(i) {eval_accuracy(stor_sample[i,], stor_labels[i,], extract_fix_param(models,i), additive = T)})
  mc_perf <- sapply(1:nrow(stor_labels), function(i) {eval_accuracy(stor_sample[i,], stor_labels[i,], get_mclust(mclust_store[[i]]))})
  
  m_m_false <- sapply(1:nrow(stor_labels), function(i) {eval_accuracy(stor_sample[i,], stor_labels[i,], extract_fix_param(models,i), additive = T, falsePos = T)})
  mc_false <-  sapply(1:nrow(stor_labels), function(i) {eval_accuracy(stor_sample[i,], stor_labels[i,], get_mclust(mclust_store[[i]]), falsePos = T)})
  
  m_m_entropy <- sapply(1:nrow(stor_labels), function(i) {eval_entropy(stor_sample[i,], stor_labels[i,], extract_fix_param(models,i), additive = T)})
  mc_entropy <-  sapply(1:nrow(stor_labels), function(i) {eval_entropy(stor_sample[i,], stor_labels[i,], get_mclust(mclust_store[[i]]))})
  
  rank_perf <- sapply(1:nrow(stor_labels), function(i){mean(rank_method(stor_sample[i,]) == stor_labels[i,])})
  rank_false <- sapply(1:nrow(stor_labels), function(i){sum(rank_method(stor_sample[i,] != stor_labels[i,] & rank_method(stor_sample[i,]) == 1))/sum(stor_labels[i,] ==0)})
  
  
  x <-
    list(
      m_m_entropy = m_m_entropy,
      mc_entropy = mc_entropy,
      m_m_perf = m_m_perf,
      mc_perf = mc_perf,
      rank_perf = rank_perf,
      m_m_false = m_m_false,
      mc_false = mc_false, 
      rank_false = rank_false,
      params = return_param_matrix,
      stan_fit = fit,
      data = stor_sample,
      labels = stor_labels
    )
  return(x)
}

sim2 <- function(mu0, mu1, theta0, sd, n=170){
  extract_fix_param <- function(float_stan, row) {
    fit_frame <- as.data.frame(float_stan)
    models <- param_vec <- apply(X = fit_frame,MARGIN = 2, mean)
    
    mu0 <- models[paste0("mu0[", row, "]")]
    mu1 <- models[paste0("mu1[", row, "]")]
    mixes <- models[paste0("theta[",row,"]")]
    sigma <- models[paste0("sigma[",row,"]")]
    
    return(list(mu0 = mu0, mu1 = mu1, sigma = sigma, mix = mixes))
  }
  
  param_matrix <- matrix(c(mu0,mu1,theta0,sd),nrow=length(mu0))
  n_mix <- sample(1:5, nrow(param_matrix), replace=T)
  
  stor_sample <- matrix(nrow=sum(n_mix), ncol=n)
  stor_labels <- matrix(nrow=sum(n_mix), ncol=n)
  current_row = 1
  return_param_matrix <- matrix(0, nrow = sum(n_mix), ncol = 4)
  
  for (i in 1:nrow(param_matrix)){
    for(m in 1:n_mix[i]) {
      temp <- generate_sample(n,mu0[i],mu1[i],sd[i],theta0[i])
      stor_sample[current_row, ] <- 
        temp[[1]]
      stor_labels[current_row,] <- 
        temp[[2]]
      return_param_matrix[current_row, ] <- param_matrix[i,]
      current_row = current_row+1
    }
  }
  
  stanFeed <- list(N = nrow(stor_sample), J = ncol(stor_sample), y = stor_sample)
  fit <- stan(file = "simple_multivar.stan", data= stanFeed, iter = 500, chains = 4 , control = list(max_treedepth = 15))
  mclust_store <- apply(stor_sample, MARGIN = 1, Mclust, G = 2, modelNames = "E")
  
  m_m_entropy <- rep(0, nrow(stor_labels))
  mc_entropy <- rep(0, nrow(stor_labels))
  
  m_m_perf <- rep(0, nrow(stor_labels))
  m_c_perf <- rep(0, nrow(stor_labels))
  m_m_perf  <- rep(0, nrow(stor_labels))
  mc_perf <- rep(0, nrow(stor_labels))
  rank_perf <- rep(0, nrow(stor_labels))
  
  param_ci_list <- as.list(rep(NA, nrow(stor_labels)))
  aij_ci_list <- as.list(rep(NA, nrow(stor_labels)))
  for(i in 1:nrow(stor_labels)) {
    print(i)
    m_m_model <- extract_fix_param(fit, i)
    mc_model <- get_mclust(mclust_store[[i]])
    
    m_m_alpha <- alpha_wrap(m_m_model, stor_sample[i,], additive = T)
    m_c_alpha <- alpha_wrap(mc_model, stor_sample[i,])
    
    m_m_entropy[i] <-  sum(-(stor_labels[i,]*log(m_m_alpha) + (1 - stor_labels[i,])*log(1 - m_m_alpha + .000001)))
    mc_entropy[i] <- sum(-(stor_labels[i,]*log(m_c_alpha) + (1 - stor_labels[i,])*log(1 - m_c_alpha + .000001)))
    
    m_m_label <- ifelse(alpha_wrap(m_m_model, stor_sample[i,], additive = T) > .5, 1, 0)
    mc_label <- ifelse(alpha_wrap(mc_model, stor_sample[i,]) > .5, 1, 0)
    
    m_m_perf[i] <- mean(m_m_label == stor_labels[i,])
    mc_perf[i] <- mean(mc_label == stor_labels[i, ])
    rank_perf[i] <- mean(rank_method(stor_sample[i,]) == stor_labels[i,])
    
    #aij_ci_list[[i]] <- lapply(stor_sample[i,] , aij_interval, stan_object = fit, row =  i)
    param_ci_list[[i]] <- param_interval(fit, i)
    
    
  }
  y <-
    list(
      m_m_entropy = m_m_entropy,
      mc_entropy = mc_entropy,
      m_m_perf = m_m_perf,
      mc_perf = mc_perf,
      rank_perf = rank_perf,
      params = return_param_matrix,
      aij_ci = aij_ci_list,
      param_ci = param_ci_list,
      stan_fit = fit,
      data = stor_sample
    )
  return(y)
}
