#' Individual testing inference for logistic regression after lasso selection
#'
#' @param X design matrix. First column is treated as an intercept.
#' @param Y.hat EM estimated response vector Y of 0s and 1s
#' @param lambda the lasso sparsity tuning parameter
#' @param beta_hat EM estimated coefficient vector
#' @param level.alpha the level of confidence
#' @return A list with the post selection confidence intervals, coefficients estimates, Truncated Normal standard deviations, and Truncated Normal limits
#'
#' @examples
#' # generate some data
#' n <- 1000
#' p <- 5
#' Se <- 0.95
#' Sp <- 0.95
#' b_true <- c(-5,3,0,1,-2,0) # first is intercept
#' X <- cbind(rep(1,n),scale(matrix(rnorm(n*p),nrow=n),TRUE,TRUE))
#'
#' pi.true <- as.numeric(logistic(X %*% b_true))
#' Y <- rbinom(n,1,logistic(pi.true))
#' Z <- rbinom(n,1,prob = Se*Y+(1-Sp)*(1-Y))
#' beta_hat <- EM_individial(X=X, Z=Z, Se=0.95, Sp=0.95, lambda=5, tol=1e-3)$beta_hat
#' Y.hat <- EM_individial(X=X, Z=Z, Se=0.95, Sp=0.95, lambda=5, tol=1e-3)$Y.hat
#' post_logreg_individual(X,Y.hat,lambda=5,beta_hat,level.alpha = 0.05)
#' @export

post_logreg_individual = function(X, Y.hat, lambda, beta_hat, level.alpha=0.05){
  pi_hat = as.numeric(logistic(X %*% beta_hat))
  W = diag(pi_hat*(1-pi_hat))
  z = X %*% beta_hat + solve(W)%*%(Y.hat-pi_hat)

  active = which(beta_hat[-1]!=0)

  k = length(active)
  s = sign(beta_hat[active+1])
  X.M = X[,c(1,(active+1))]

  inv_info = solve(t(X.M) %*% W %*% X.M)
  A1 = -diag(s)
  b1 = -diag(s) %*% inv_info[-1,] %*% c(0, lambda*s)

  theta_bar = inv_info %*% t(X.M) %*% W %*% z
  alpha_bar = theta_bar[1]; beta_bar = theta_bar[-1]

  H = -t(X.M) %*% W %*% X.M + t(X.M) %*% diag(Y.hat*(1-Y.hat)) %*% X.M #Hessian for individual testing
  inv_H = solve(-H)
  Sigma = inv_H[-1,-1]

  CI = matrix(0,k,2)
  pvalue = c(); p_value = c()
  TNsd = c()
  limits = matrix(0,k,2)
  for (j in 1:k){
    ej = matrix(0,k,1); ej[j] = 1
    eta = ej
    c = Sigma %*% eta %*% solve(t(eta)%*%Sigma%*%eta)
    r = beta_bar - c %*% t(eta) %*% beta_bar
    A1r = A1 %*% r
    A1c = A1 %*% c

    quantity = c()
    for (i in 1:length(A1c)) {
      quantity[i] = (b1[i]-A1r[i]) / as.vector(A1c[i])
    }
    V.neg = suppressWarnings(max(quantity[A1c<0]))
    V.pos = suppressWarnings(min(quantity[A1c>0]))

    sd =  sqrt(Sigma[j,j])
    CI.L = function(x){ ptruncnorm(as.vector(t(eta) %*% beta_bar), a=V.neg, b=V.pos, mean=x, sd=sd) + level.alpha/2-1 }
    CI.U = function(x){ ptruncnorm(as.vector(t(eta) %*% beta_bar), a=V.neg, b=V.pos, mean=x, sd=sd) - level.alpha/2 }
    grid = seq(-10,10,0.001)
    CI.LOW = grid[which.min(abs(CI.L(grid)))]
    CI.UP = grid[which.min(abs(CI.U(grid)))]
    CI[j,] = round(c(CI.LOW, CI.UP),4)
    limits[j,] = round(c(V.neg,V.pos),4)
    TNsd[j] = round(TG.limits(beta_bar, A1, b1, eta, Sigma)$sd,4)
  }
  return(list(CI=CI,beta_bar=beta_bar,TNsd=TNsd,limits=limits))
}


#' Group testing inference for logistic regression after lasso selection
#'
#' @param X design matrix. First column is treated as an intercept.
#' @param Z observed group response vector of 0s and 1s
#' @param Y.hat EM estimated response vector of 0s and 1s
#' @param Se sensitivity
#' @param Sp specificity
#' @param lambda the lasso sparsity tuning parameter
#' @param beta_hat EM estimated coefficient vector
#' @param n_groups number of groups
#' @param level.alpha the level of confidence
#' @return A list with the post selection confidence intervals, coefficients estimates, Truncated Normal standard deviations, and Truncated Normal limits
#'
#' @examples
#' # generate some data
#' n <- 1000
#' p <- 5
#' Se <- 0.95
#' Sp <- 0.95
#' b_true <- c(-5,3,0,1,-2,0) # first is intercept
#' X <- cbind(rep(1,n),scale(matrix(rnorm(n*p),nrow=n),TRUE,TRUE))
#' pi.true <- as.numeric(logistic(X %*% b_true))
#' Y <- rbinom(n,1,prob = pi.true)
#'
#' n_groups <- 250
#' m <- n/n_groups #number of samples within a group
#' Z <- c()
#' for(j in 1:n_groups){
#'   maxY <- max(Y[((j-1)*m+1):(j*m)])
#'   Z[j] <- rbinom(1,1,prob=Se*maxY+(1-Sp)*(1-maxY))
#' }
#'
#' EM_result <- EM_group(X=X, Z=Z, Se=0.95, Sp=0.95, lambda=5, tol=1e-3, n_groups=250)
#' beta_hat <- EM_result$beta_hat
#' Y.hat <- EM_result$Y.hat
#' post_logreg_group(X,Z,Y.hat,Se,Sp,lambda=5,beta_hat,n_groups,level.alpha = 0.05)
#' @export

post_logreg_group = function(X, Z, Y.hat, Se, Sp, lambda, beta_hat, n_groups, level.alpha=0.05){
  n = nrow(X)
  p = ncol(X) - 1
  m = n/n_groups
  pi_hat = as.numeric(logistic(X %*% beta_hat))
  W = diag(pi_hat*(1-pi_hat))
  z = X %*% beta_hat + solve(W)%*%(Y.hat-pi_hat)

  active = which(beta_hat[-1]!=0)

  n_active = length(active)
  s = sign(beta_hat[active+1])
  X.M = X[,c(1,(active+1))]

  inv_info = solve(t(X.M) %*% W %*% X.M)
  A1 = -diag(s)
  b1 = -diag(s) %*% inv_info[-1,] %*% c(0, lambda*s)

  theta_bar = inv_info %*% t(X.M) %*% W %*% z
  beta_bar = theta_bar[-1]

  # obtain Hessian
  XcovX = matrix(0,n_active+1,n_active+1)
  for(j in 1:n_groups){
    X_Pj = X.M[((j-1)*m+1):(j*m),]
    pi_j = pi_hat[((j-1)*m+1):(j*m)]
    XcovX = XcovX + t(X_Pj) %*% cov_matrix(Se = Se, Sp = Sp, pi = pi_j, Z=Z[j]) %*% X_Pj
  }

  H = -t(X.M) %*% W %*% X.M + XcovX

  inv_H = solve(-H)
  Sigma = inv_H[-1,-1]

  CI = matrix(0,n_active,2)
  TNsd = c()
  limits = matrix(0,n_active,2)
  for (j in 1:n_active){
    ej = matrix(0,n_active,1); ej[j] = 1
    eta = ej
    c = Sigma %*% eta %*% solve(t(eta)%*%Sigma%*%eta)
    r = beta_bar - c %*% t(eta) %*% beta_bar
    A1r = A1 %*% r
    A1c = A1 %*% c

    quantity = c()
    for (i in 1:length(A1c)) {
      quantity[i] = (b1[i]-A1r[i]) / as.vector(A1c[i])
    }
    V.neg = suppressWarnings(max(quantity[A1c<0]))
    V.pos = suppressWarnings(min(quantity[A1c>0]))
    V.0 = suppressWarnings(min((b1-A1r)[A1c==0]))

    sd = sqrt(Sigma[j,j])
    CI.L = function(x){ ptruncnorm(as.vector(t(eta) %*% beta_bar), a=V.neg, b=V.pos, mean=x, sd=sd) + level.alpha/2-1 }
    CI.U = function(x){ ptruncnorm(as.vector(t(eta) %*% beta_bar), a=V.neg, b=V.pos, mean=x, sd=sd) - level.alpha/2 }
    grid = seq(-10,10,0.001)
    CI.LOW = grid[which.min(abs(CI.L(grid)))]
    CI.UP = grid[which.min(abs(CI.U(grid)))]
    CI[j,] = round(c(CI.LOW, CI.UP),4)
    limits[j,] = round(c(V.neg,V.pos),4)
    TNsd[j] = round(TG.limits(beta_bar, A1, b1, eta, Sigma)$sd,4)
  }
  return(list(CI=CI,beta_bar=beta_bar,TNsd=TNsd,limits=limits))
}




#' Classical inference for logistic regression after lasso selection, individual testing
#'
#' @param X.M design matrix with selected columns. First column is treated as an intercept.
#' @param Y.hat EM estimated response vector of 0s and 1s
#' @param beta_hat EM estimated coefficient vector
#' @param level.alpha the level of confidence
#' @return A list with the Naive inference confidence intervals and p values
#'
#' @export

naive_logreg_individual = function(X.M, Y.hat, beta_hat, level.alpha=0.05){
  pi_hat = as.numeric(logistic(X.M %*% beta_hat))
  W = diag(pi_hat*(1-pi_hat))

  k = length(beta_hat)-1

  H = -t(X.M) %*% W %*% X.M + t(X.M) %*% diag(Y.hat*(1-Y.hat)) %*% X.M #Hessian for individual testing
  inv_H = solve(-H)
  Sigma = inv_H[-1,-1]

  CI = matrix(0,k,2)
  p_value = c()
  for (j in 1:k){
    se = sqrt(Sigma[j,j])
    test_stats = beta_hat[j+1]/se
    p_value[j] = 2*pnorm(-abs(test_stats),lower.tail = TRUE)
    CI.LOW = beta_hat[j+1] - qnorm(level.alpha/2,lower.tail=FALSE)*se
    CI.UP = beta_hat[j+1] + qnorm(level.alpha/2,lower.tail=FALSE)*se
    CI[j,] = round(c(CI.LOW, CI.UP),4)
  }

  return(list(CI=CI,pvalue=round(p_value,4)))
}


#' Classical inference for logistic regression after lasso selection, group testing
#'
#' @param X.M design matrix with selected columns. First column is treated as an intercept.
#' @param Z observed group response vector of 0s and 1s
#' @param Se sensitivity
#' @param Sp specificity
#' @param beta_hat EM estimated coefficient vector
#' @param level.alpha the level of confidence
#' @param n_groups number of groups
#' @return A list with the Naive inference confidence intervals and p values
#'
#' @export

naive_logreg_group = function(X.M, Z, Se, Sp, beta_hat, level.alpha=0.05, n_groups){
  n = nrow(X.M)
  m = n/n_groups
  n_active = ncol(X.M)-1
  pi_hat = as.numeric(logistic(X.M %*% beta_hat))
  W = diag(pi_hat*(1-pi_hat))
  k = length(beta_hat)-1

  # obtain Hessian
  XcovX = matrix(0,n_active+1,n_active+1)
  for(j in 1:n_groups){
    X_Pj = X.M[((j-1)*m+1):(j*m),]
    pi_j = pi_hat[((j-1)*m+1):(j*m)]
    XcovX = XcovX + t(X_Pj) %*% cov_matrix(Se = Se, Sp = Sp, pi = pi_j, Z=Z[j]) %*% X_Pj
  }

  H = -t(X.M) %*% W %*% X.M + XcovX

  inv_H = solve(-H)
  Sigma = inv_H[-1,-1]

  CI = matrix(0,k,2)
  p_value = c()
  for (j in 1:k){
    se = sqrt(Sigma[j,j])
    test_stats = beta_hat[j+1]/se
    p_value[j] = 2*pnorm(-abs(test_stats),lower.tail = TRUE)
    CI.LOW = beta_hat[j+1] - qnorm(level.alpha/2,lower.tail=FALSE)*se
    CI.UP = beta_hat[j+1] + qnorm(level.alpha/2,lower.tail=FALSE)*se
    CI[j,] = round(c(CI.LOW, CI.UP),4)
  }

  return(list(CI=CI,pvalue=round(p_value,4)))
}
