logistic = function(x){1/(1+exp(-x))}


#' Function to compute logistic regression parameters with a lasso penalty using IRLS
#'
#' @param X design matrix. First column is treated as an intercept.
#' @param Y response vector of 0s and 1s
#' @param lambda the lasso sparsity tuning parameter
#' @param tol convergence tolerance
#'
#' @return numeric vector of estimated regression coefficients
#'
#' @examples
#' # generate some data
#' n <- 1000
#' p <- 5
#' b_true <- c(0,3,-2,1,0,0) # first is intercept
#' X <- cbind(rep(1,n),matrix(rnorm(n*p),nrow=n))
#' eta <- X %*% b_true
#' Y <- rbinom(n,1,logistic(eta))
#'
#' logreg_lasso(X,Y,lambda = 5)
#' @export

logreg_lasso <- function(X,Y,lambda,tol=1e-5){
  p <- ncol(X)
  b <- rep(0,p)
  conv0 <- F
  while(conv0 == F){
    b00 <- b
    eta <- X %*% b
    pr <- logistic(eta)
    w <- pr*(1-pr)
    z <- eta + (Y - pr)/w
    conv1 <- F
    while(conv1 == F){
      b0 <- b
      # update intercept
      b[1] <- sum(w*(z - X[,-1] %*% b[-1]))/sum(w)
      # cycle through the rest of the coefficients
      for(j in 2:p){
        zj <- z - X[,-j] %*% b[-j] # there is a trick to make this faster...
        bj <- sum(X[,j]*w*zj)
        # soft-thresholding
        if(bj > lambda){
          bj <- bj - lambda
        } else if(bj < -lambda){
          bj <- bj + lambda
        } else {
          bj <- 0
        }
        b[j] <- bj / sum(X[,j]^2 * w)
      }
      conv1 <- max(abs(b - b0)) < tol
    }
    conv0 <- max(abs(b - b00)) < tol
  }
  output <- list(b = b,
                 lambda = lambda)
  return(output)
}
