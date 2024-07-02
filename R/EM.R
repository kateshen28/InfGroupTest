#' Function to perform EM algorithm for individual testing
#'
#' @param X design matrix. First column is treated as an intercept.
#' @param Z observed individual response vector of 0s and 1s
#' @param Se sensitivity
#' @param Sp specificity
#' @param lambda the lasso sparsity tuning parameter
#' @param tol convergence tolerance
#'
#' @return numeric vector of estimated regression coefficients beta_hat, and estimated response Y.hat
#'
#' @examples
#' # generate some data
#' n <- 1000
#' p <- 5
#' Se <- 0.95
#' Sp <- 0.95
#' beta_true <- c(0,3,-2,1,0,0) # first is intercept
#' X <- cbind(rep(1,n),matrix(rnorm(n*p),nrow=n))
#' pi.true <- as.numeric(logistic(X %*% beta_true))
#' Y <- rbinom(n,1,prob = pi.true)
#' Z <- rbinom(n,1,prob = Se*Y+(1-Sp)*(1-Y))
#'
#' EM_individial(X=X, Z=Z, Se=0.95, Sp=0.95, lambda=5, tol=1e-3)
#' @export

EM_individial = function(X, Z, Se, Sp, lambda, tol=1e-3){
  n = nrow(X); p = ncol(X)-1
  beta_hat = rep(1,p+1)
  Y.hat = rep(0,n)
  error = 1
  while(error>tol){
    pi = as.numeric(logistic(X %*% beta_hat))
    Y.hat_old = Y.hat
    Y.hat = ifelse(Z == 1,Se*pi/(Se*pi + (1-Sp)*(1-pi)),(1-Se)*pi/((1-Se)*pi + Sp*(1-pi)))

    # coordinate descent
    beta_hat = logreg_lasso(X,Y.hat,lambda)$b

    #tol = max( abs(Y.hat - Y.hat_old) )
    error = sum((Y.hat - Y.hat_old)^2)
  }
  return(list(beta_hat=beta_hat, Y.hat=Y.hat))
}


#' Function to calculate covariance matrix given Z for pool Pj
#'
#' @param Se sensitivity
#' @param Sp specificity
#' @param pi a vector of probabilities of length mj, number of samples within jth pool
#' @param Z observed group response either 0 or 1
#'
#' @return a square matrix of size mj*mj where mj is number of samples within jth pool
#'
#' @examples
#' # generate some data
#' n <- 1000
#' p <- 5
#' n_groups <- 250
#' m <- n/n_groups
#' Se <- 0.95
#' Sp <- 0.95
#' X <- cbind(rep(1,n),matrix(rnorm(n*p),nrow=n))
#' beta_hat <- c(-5,3,-2,1,0,0) # first is intercept
#' pi.hat <- as.numeric(logistic(X %*% beta_hat))
#'
#' cov_matrix(Se=0.95, Sp=0.95, pi=pi.hat[1:m], Z=0) #first pool
#' @export

cov_matrix = function(Se, Sp, pi, Z){
  mj = length(pi)
  p_z0 = (1-Se)-(1-Se)*prod(1-pi)+Sp*prod(1-pi)
  p_z1 = Se-Se*prod(1-pi)+(1-Sp)*prod(1-pi)
  cov = matrix(0,mj,mj)
  for (i in 1:(mj-1) ){
    for (k in (i+1):mj){
      if (Z == 0){
        cov[i,k] = ( (1-Se)*pi[i]*pi[k]*(p_z0-1+Se) ) / p_z0^2
        cov[k,i] = cov[i,k]
      }
      if (Z == 1){
        cov[i,k] = ( Se*pi[i]*pi[k]*(p_z1-Se) ) / p_z1^2
        cov[k,i] = cov[i,k]
      }
    }
  }
  for(i in 1:mj){
    if (Z == 0){cov[i,i] = (1-Se)*pi[i]/p_z0-(1-Se)^2*pi[i]^2/p_z0^2}
    if (Z == 1){cov[i,i] = Se*pi[i]/p_z1-Se^2*pi[i]^2/p_z1^2}
  }
  return(cov)
}


E.Y = function(Se, Sp, pi, Z){
  mj = length(pi)
  res = rep(0,mj)
  if (Z == 0){res = (1-Se)*pi/( (1-Se)-(1-Se)*prod(1-pi)+Sp*prod(1-pi) )}
  if (Z == 1){res = Se*pi/( Se-Se*prod(1-pi)+(1-Sp)*prod(1-pi) )}
  return(res)
}

#' Function to perform EM algorithm for group testing
#'
#' @param X design matrix. First column is treated as an intercept.
#' @param Z observed group response vector of 0s and 1s
#' @param Se sensitivity
#' @param Sp specificity
#' @param lambda the lasso sparsity tuning parameter
#' @param tol convergence tolerance
#'
#' @return numeric vector of estimated regression coefficients beta_hat, and estimated response Y.hat
#'
#' @examples
#' # generate some data
#' n <- 1000
#' p <- 5
#' Se <- 0.95
#' Sp <- 0.95
#' beta_true <- c(-5,3,0,1,-2,0) # first is intercept
#' X <- cbind(rep(1,n),matrix(rnorm(n*p),nrow=n))
#' pi.true <- as.numeric(logistic(X %*% beta_true))
#' Y <- rbinom(n,1,prob = pi.true)
#' n_groups <- 250
#' m <- n/n_groups
#' Z <- c()
#' for (j in 1:n_groups){
#'   maxY <- max(Y[((j-1)*m+1):(j*m)])
#'   Z[j] <- rbinom(1,1,prob = Se*maxY+(1-Sp)*(1-maxY))
#' }
#'
#' EM_group(X=X, Z=Z, Se=0.95, Sp=0.95, lambda=5, tol=1e-3, n_groups=250)
#' @export

EM_group = function(X, Z, Se, Sp, lambda, tol=1e-3, n_groups){
  n = nrow(X); p = ncol(X)-1
  m = n/n_groups #number of samples within a group
  beta_hat = rep(1,p+1)
  Y.hat = rep(0,n)
  error = 1
  while(error>tol){
    pi = as.numeric(logistic(X %*% beta_hat))
    Y.hat_old = Y.hat
    for(j in 1:n_groups){
      Y.hat[((j-1)*m+1):(j*m)] = E.Y(Se,Sp,pi[((j-1)*m+1):(j*m)],Z[j])
    }

    # coordinate descent
    beta_hat = logreg_lasso(X,Y.hat,lambda)$b

    #tol = max( abs(Y.hat - Y.hat_old) )
    error = sum((Y.hat - Y.hat_old)^2)
  }
  return(list(beta_hat=beta_hat, Y.hat=Y.hat))
}

