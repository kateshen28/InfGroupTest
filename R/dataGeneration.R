#' Function to generate individual testing data
#'
#' @param n number of observations of X
#' @param p number of variables of X
#' @param Se sensitivity
#' @param Sp specificity
#' @param beta_true true regression coefficient
#' @return A list of covariates X and observed responses Z
#'
#' @examples
#' # generate individual testing data
#' n <- 1000
#' p <- 5
#' Se <- 0.95
#' Sp <- 0.95
#' beta_true <- c(0,3,-2,1,0,0) # first is intercept
#' generate.data(n, p, Se, Sp, beta_true)
#' @export
generate.data = function(n, p, Se, Sp, beta_true){
  X = cbind(rep(1,n),scale(matrix(rnorm(n*p),nrow=n),T,T))
  pi.true = as.numeric(logistic(X %*% beta_true))
  Y = rbinom(n,1,prob = pi.true)
  Z = rbinom(n,1,prob=Se*Y+(1-Sp)*(1-Y))
  output=list(X=X,Z=Z)
  return(output)
}

#' Function to generate group testing data
#'
#' @param n number of observations of X
#' @param p number of variables of X
#' @param Se sensitivity
#' @param Sp specificity
#' @param beta_true true regression coefficient
#' @param n_groups number of groups
#' @return A list of covariates X and observed responses Z
#'
#' @examples
#' # generate group testing data
#' n <- 1000
#' p <- 5
#' Se <- 0.95
#' Sp <- 0.95
#' beta_true <- c(0,3,-2,1,0,0) # first is intercept
#' n_groups <- 250
#' generate.data.group(n, p, Se, Sp, beta_true, n_groups)
#' @export
generate.data.group = function(n, p, Se, Sp, beta_true, n_groups){
  X = cbind(rep(1,n),scale(matrix(rnorm(n*p),nrow=n),T,T))
  pi.true = as.numeric(logistic(X %*% beta_true))
  Y = rbinom(n,1,prob = pi.true)

  m = n/n_groups #number of samples within a group
  Z = c()
  for(j in 1:n_groups){
    maxY = max(Y[((j-1)*m+1):(j*m)])
    Z[j] = rbinom(1,1,prob=Se*maxY+(1-Sp)*(1-maxY))
  }
  output=list(X=X,Z=Z)
  return(output)
}
