#' example2: Create the data used in the second setting of the simulation section of Li, Goldberg, \& Zheng
#'
#' See details in the simulation section
#'
#' @param seed  The seed for the simulation.
#' @param n The number of observations.
#' @param var.est variance of the score.
#' @param thres risk threshold for calculating the age to start screening
#' @param upper.age  the upper time denoted by t in the paper
#' @param censoring   Either "simple" or "complex". Simple means that the censoring time is independent of the failure time. Complex means that it dependent.
#'
#' @return A data frame with covariate 'Z',
#' observed time 'time'=Z,
#' the cometing risk 'status' which is 0 for censoring, 1 for the event of interest, and 2 for the competing risk=status,
#' the censoring time 'C', event of interest time 'T1', competing risk time 'T2',
#' and 'tau1' and 'tau2' which are the lower and upper times denoted by S and t in the paper
#'
#' @import tibble
#' @importFrom stats integrate dnorm rnorm rexp pnorm qnorm runif uniroot approxfun
#' @export
example2 <-  function(seed = 1, n = 1000, var.est = 0.1, thres=0.001, upper.age = 60, censoring="simple"){
# flag: TRUE (censoring not depending on Z) False otherwise.
# var.est: variance of the score
# thres: risk threshold for calculating the age to start screening. If it is equal to 0, starting
#        age is set to be 50, corresponding to the recommended age to start screening
  if (censoring == "simple") {flag<-TRUE} else {flag<-FALSE}
  dat.info <-  create.data.crc(n=n, var.est=var.est, seed=seed, cen.flag = flag)
  tau1 <-  calc.risk(Delta=10, dat.info,thres)
  dat <-  cbind(dat.info$dat,tibble(tau1=tau1,tau2=upper.age)) %>% tibble::as_tibble()
  return(dat)
}

create.data.crc <-  function(n, var.est, seed, cen.flag = TRUE){
#########################################################################
#  composit CRC incidence and mortality from other causes               #
#  CRC incidence in SEER13 (1992-2005)                                  #
#  var.est: variance of the score                                       #
#  cen.flag=TRUE, censoring does not depend on Z and FALSE otherwise    #
#########################################################################
  set.seed(seed)

#### CRC incidence rates, mortality rates and current age distribution


  lambda<- CRCMortilityAge

#### CRC and finding the baseline hazard function for CRC

  surv.crc = exp(-cumsum(lambda$CRC.HzaRate))   ## Survival probability (x2, to see more common disease)

  find.basesurv <-  function(a, surv.t, var.est){
    ##########################
    ### a: baseline survival probability
    ### surv.t: survival probability at time t
    ### var.est: variance of risk score. Note the mean of risk score is assumed to be 0
    ##########################
    fcn = function(x, a, var.est){
      return(a^(exp(x))*stats::dnorm(x,mean = 0, sd = var.est^0.5))
    }
    return(surv.t - stats::integrate(fcn, lower = -Inf, upper = Inf, a=a, var.est=var.est)$value)
  }

  surv0.crc = vector()   ### S(t|Z=0)
  for (i in 1:length(surv.crc)){
    if (surv.crc[i]==1){surv0.crc[i]=1} else{
      surv0.crc[i] <- stats::uniroot(find.basesurv, c(0.0001,1), tol = 0.000001, surv.t = surv.crc[i], var.est)$root
    }
  }

  lambda0 = c(0, diff(-log(surv0.crc))) # Baseline hazard function for CRC

#### Risk score
  Z = stats::rnorm(n,mean=0,sd= var.est^0.5)   ### risk score

#### Age at disease diagnosis

  u1 = stats::runif(n)^exp(-Z)
  S.crc <- stats::approxfun(surv0.crc, 1:110, method="linear", rule=2)
  T1 = S.crc(u1)

#### Competing risks of death

  surv.mort <- exp(-cumsum(lambda$Mortality.HzdRate))
  S.mort <- stats::approxfun(surv.mort, 1:110, method="linear", rule=2)
#  T2 = S.mort(runif(n)^(exp(-log(2)*Z)))
  T2 = S.mort(runif(n)^(exp(-Z)))

#####  Censoring from the current age

  S3 = stats::approxfun(lambda$CumProb.CurrentAge, 1:110, method="linear", rule=2)
  if (cen.flag) { age = S3(runif(n)) } else {
    age = S3(runif(n, 0, 1)^exp(-Z))  ### censoring age depends on Z
  }


#### Combining the event of interest (T1.age), of competing risks (T2.age), and censoring (C)
  T1 = round(T1, digits = 0)
  T2 = round(T2, digits = 0)
  age = round(age, digits = 0)
  X = pmin(T1, T2, age)

#### if status = 0, censoring, 1, event of interest, 2 event of competing risks
  status = rep(0, n)
  status[X==T2] = 2
  status[X==T1] = 1

  dat = tibble::tibble(T1 = T1, T2 = T2, C=age, time=X, status=status, Z = Z)
  dat.info = list(dat = dat, lambda0T1 = lambda0, lambda0T2 = lambda$Mortality.HzdRate)
  return(dat.info)
}


risk.estimation = function(age0, Delta, dat.info){
#########################################################################
#### lambda0T1: baseline hazard function for the event of interest
#### lambda0T2: (baseline) hazard for the competing risks event
#### score: exp(Z) for T1
#### score2: exp(Z) for T2
#### Calculation the Delta year risk at any given age0
##########################################################################

  n = nrow(dat.info$dat)
  lambda0T1 = dat.info$lambda0T1
  lambda0T2 = dat.info$lambda0T2
  score = exp(dat.info$dat$Z)
  if ("Z2" %in% names(dat.info$dat)) {score2 = exp(dat.info$dat$Z2)} else {score2 = rep(1,n)}
  Lambda0 = cumsum(c(0,lambda0T1[(age0+1):(age0+Delta-1)]))
  cumlambda0 = exp(-Lambda0%*%t(score) - cumsum(c(0,lambda0T2[(age0+1):(age0+Delta-1)]))%*%t(score2)) ### with competing risks
  risk = apply((lambda0T1[(age0+1):(age0+Delta)]%*%t(score))*cumlambda0, 2, sum)
  cumlambda0 = exp(-Lambda0%*%t(score))    ### without competing risks
  risk.nocompeting = apply((lambda0T1[(age0+1):(age0+Delta)]%*%t(score))*cumlambda0, 2, sum)
  return(tibble(risk=risk, risk.nocompeting = risk.nocompeting))
}

max.age = function(risk, thres){
##############################################################################
### calculating the minimum age at which the risk exceeds the threshold
### if thres==0, set minimum age = 50. In the CRC example, 50 corresponds to the standard
###    starting age for screening
##############################################################################

  if (thres>0) {
    temp = which(risk>thres)
    if (length(temp)==0) {ss=length(risk)} else {ss=min(temp)}
  }
  if (thres==0) {ss = 50}
  return(ss)
}

calc.risk <- function(Delta,dat.info,thres){
# Delta: Delta year of risk
  n <- nrow(dat.info$dat)
  risk = matrix(0,n,100)   ## with competing risk
  risk0 = matrix(0,n,100)  ## with no competing risk
  for (j in 1:100){
    temp = risk.estimation(j,Delta, dat.info)
    risk[,j] = temp$risk
    risk0[,j] = temp$risk.nocompeting
  }
  start.age = apply(risk,1,max.age,thres=thres)
  return(start.age)
}

