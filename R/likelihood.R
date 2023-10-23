#' Log-likelihood function for logistic regression after lasso selection, individual testing
#'
#' @param Z observed group response vector of 0s and 1s
#' @param Se sensitivity
#' @param Sp specificity
#' @param pi_hat a vector of estimated probabilities of length mj
#' @return log likelihood
#'
#' @export
loglike_ind = function(Z, Se, Sp, pi_hat){
  l = sum( Z*log( Se*pi_hat+(1-Sp)*(1-pi_hat) ) + (1-Z)*log( Sp*(1-pi_hat)+(1-Se)*pi_hat ) )
}


#' Log-likelihood function for logistic regression after lasso selection, group testing
#'
#' @param Z observed group response vector of 0s and 1s
#' @param Se sensitivity
#' @param Sp specificity
#' @param pi_hat a vector of estimated probabilities of length mj
#' @param n_groups number of groups
#' @return log likelihood
#'
#' @export
loglike_group = function(Z, Se, Sp, pi_hat, n_groups, m){
  l = 0
  for(j in 1:n_groups){
    pi_j = pi_hat[((j-1)*m+1):(j*m)]
    part1 = Z[j]*log( Se - Se*prod(1-pi_j) + (1-Sp)*prod(1-pi_j))
    part2 = (1-Z[j])*log( (1-Se) - (1-Se)*prod(1-pi_j) + Sp*prod(1-pi_j) )
    l = l + part1 + part2
  }
}
