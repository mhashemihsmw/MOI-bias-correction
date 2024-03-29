#' Derives the probability of pathological data
#'
#' @param N integer; sample size 
#' @param lambda numeric; MOI parameter
#' @param p numeric vector; lineage frequencies
#'
#' @return numeric; probability of pathological data
#' 
#' @export
#'
#' @examples
#' \donotrun{
#' prob_pathological(150, 1.23, c(0.64, 0.36))
#' }
prob_pathological <- function(N, lambda, p){
    lep <- lambda*p
    dk <- exp(lep) - 1
    d <- sum(dk)
    d0 <- 1/(exp(lambda) - 1)
    x <- (1 - d*d0) 
    y <- N*(d0 + 1)
    den <- y*x
    nom<- (d0 + 1/2)*d - d0*((d^2) - sum(dk^2))/(2*x)
    
    q1 <- (d*d0)^N 
    q2 <- sum((dk*d0)^N)
    q3 <- (1 - prod(1 - (1 - exp(-lep))^N))/(1 - exp(-lambda))^N
    
    q1[is.nan(q1)==T] <- 0
    q2[is.nan(q2)==T] <- 0
    q3[is.nan(q3)==T] <- 0
    
    p_pathologic <- q3 + q1 - q2
    p_pathologic
}
