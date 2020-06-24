#' example1: Create the data used in the first setting of the simulation section of Li, Goldberg, \& Zheng
#'
#' See details in the simulation section
#'
#' @param seed  The seed for the simulation
#' @param n The number of observations
#' @param censoring   Either "simple" or "complex". Simple means that the censoring time is independent of the failure time. Complex means that it dependent.
#' @param to_round Boolean, If true, times are rounded down to two decimal places.
#'
#' @return A data framewith covariate 'Z',
#' observed time 'time'=Z,
#' the cometing risk 'status' which is 0 for censoring, 1 for the event of interest, and 2 for the competing risk=status,
#' the censoring time 'C', event of interest time 'T1', competing risk time 'T2',
#' and 'tau1' and 'tau2' which are the lower and upper times denoted by S and t in the paper
#'
#' @import tibble
#' @importFrom stats rnorm rexp pnorm qnorm runif uniroot
#' @export
example1 <- function(seed=1,n=500,censoring = "simple",
                     to_round=T){
  set.seed(seed = seed)


  lambda1 <- 1.8
  lambda2 <- 1.4
  lambda0 <- 3

  Z <- stats::rnorm(n=n,sd = 1)
  tau1 <- floor_dec(exp(-lambda1*exp(Z)),2)
  tau2 <- 0.3

  T1 <- stats::rexp(n,rate = exp(Z)*lambda1)
  T2 <- stats::rexp(n,rate = exp(0.5*Z)*lambda2)
  if(censoring=="simple"){
    C <- stats::rexp(n,rate = lambda0)+0.01
  } else {
    C<- stats::rexp(n,rate = exp(2*Z)*lambda0)+0.01
  }


  if(to_round){
    T1 <-  floor_dec(T1,2)
    T2 <- floor_dec(T2,2)
    C <- floor_dec(C,2)
  }
  time <-  pmin(C,T1,T2)
  status <- ifelse(time==T1,1,ifelse(time==T2,2,0))

  df <- tibble::tibble(Z=Z,time= time,
               status=status,
               C=C,T1=T1,T2=T2,tau1 = tau1,tau2 = tau2, censoring = censoring)
  return(df)
}

