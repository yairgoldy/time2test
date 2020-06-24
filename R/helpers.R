

floor_dec <- function(x, level=1) round(x - 5*10^(-level-1), level)
ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)



#' @importFrom arules discretize
#' @importFrom stringr str_locate str_sub str_length
discretize_Z_fun <- function(Z,breaks=40){
  if(length(unique(Z))<30)  return(Z)


  Zdist <- arules::discretize(Z,breaks = breaks)
  lZ <- levels(Zdist)
  comma_place <- stringr::str_locate(lZ,",")[,1]
  lower <-  stringr::str_sub(lZ,2,comma_place-1) %>% as.numeric()
  upper <-  stringr::str_sub(lZ,comma_place+1,stringr::str_length(lZ)-1) %>% as.numeric()
  means <- (lower+upper)/2
  ZZ <- dplyr::left_join(tibble(Z=Z,Zdist=as.character(Zdist)), tibble::tibble(Zdist=lZ,Zmean=means),by = "Zdist")
  return(ZZ$Zmean)
}

#' @import dplyr
#' @importFrom survival survfit Surv
#' @importFrom stats approxfun
KM_Censoring <- function(df,eps=1e-12){


  KM <- summary(survival::survfit(survival::Surv(time,status==0) ~ 1,  type="kaplan-meier",  data=df))

    # The kaplan meier estimator of the censoring
  # shifted so that at a point time T KMC(T)=G(T-)
  KM <- tibble::tibble(time=KM$time,surv=KM$surv,n.event=KM$n.event,n.risk=KM$n.risk,CN=cumsum(n.event)) %>%
    dplyr::filter(surv!=0)

  KMC <- stats::approxfun( KM$time+eps, KM$surv,method="constant",
                    yleft = 1,yright = KM$surv[length(KM$surv)])
  # The at risk process
  R <- stats::approxfun( KM$time-eps, KM$n.risk,method="constant",
                  yleft = KM$n.risk[1],yright = KM$n.risk[length(KM$n.risk)])

  # dN process
  N <- stats::approxfun( KM$time,KM$CN ,method="constant",
                  yleft = 0,yright =KM$CN[length(KM$CN)] )
  KM <- dplyr::mutate(KM,LN = cumsum(n.event/n.risk))
  Lambda <- stats::approxfun( KM$time,KM$LN ,method="constant",
                       yleft = 0,yright =KM$LN[length(KM$LN)] )
  Censoring_times <- KM$time



  return( list(KMC=KMC,R=R,N=N,Lambda=Lambda,Censoring_times=Censoring_times))
}


CIF <- function(Z,competing_risk_model,event=1){

  if(event == 1){
    crm <- dplyr::mutate(competing_risk_model,
                  S=exp(-exp(Z*b1)*cumhaz1-exp(Z*b2)*cumhaz2),
                  lam1=exp(Z*b1)*(diff(c(0,cumhaz1))),
                  lam1prodS=lam1*S,
                  CIF=cumsum(lam1prodS))
  }
  if (event ==2) {
    crm <- dplyr::mutate(competing_risk_model,
                  S=exp(-exp(Z*b1)*cumhaz1-exp(Z*b2)*cumhaz2),
                  lam2=exp(Z*b2)*(diff(c(0,cumhaz2))),
                  lam2prodS=lam2*S,
                  CIF=cumsum(lam2prodS))

  }

  if (event==3){
    crm <-dplyr:: mutate(competing_risk_model,
                  CIF=1-exp(-exp(Z*b1)*cumhaz1-exp(Z*b2)*cumhaz2))


  }
  return(stats::approxfun( crm$time,crm$CIF, method="constant",
                    yleft = 0,yright = crm$CIF[length(crm$CIF)]))
}

SgivenZ <- function(Znew,competing_risk_model){
  return( function(x){1-CIF(Znew,competing_risk_model,event = 3)(x)})
}




sumQZiCk_RiCk <- function(Ct,df,Qvals){
  res <- dplyr::mutate(df,RiCk= as.numeric(time>=Ct))
  QvalCk <- dplyr::filter(Qvals,Ck==Ct) %>%
    dplyr::select(discZ,QZiCk)
  res <- dplyr::left_join(res,QvalCk,by="discZ") %>%
    dplyr::mutate(QZiCk_RiCk=QZiCk*RiCk) %>%
    dplyr::summarise(sumQR=sum(QZiCk_RiCk)) %>%
    unlist(use.names = F)
  return(res)
}






calc_event <- function(df,type=1){

  switch (as.character(type),
          "1" = {return( df$tau1 < df$time & df$time <= df$tau2 &
                           df$status == 1)
          },
          "2" = {return( df$tau1 < df$time & df$time <= df$tau2 &
                           df$status == 2)
          },
          "3" = {return( -( df$time <= df$tau2  & df$tau1<= df$tau2 &
                              df$status !=0))
          },
          "4" = {return( (df$time<df$tau1) & df$time <= df$tau2 &
                           df$status == 1)
          },
          "5" = {return( (df$time<df$tau1) & df$time <= df$tau2 &
                           df$status == 2)
          },
          "6" = {return(- ( df$time <= df$tau2  & df$tau1> df$tau2 &
                              df$status !=0))
          }

  )
}

calc_event_t <- function(df,t, type = 1){
  return( ( calc_event(df, type) & df$time>=t ) )
}

calc_W <- function(Ck,df, type = 1){
  mt = calc_event_t(df,Ck,type)
  Wk <- dplyr::mutate(df,mt=mt,val=mt/G) %>%
    dplyr::summarise(res=sum(val)) %>%
    unlist(use.names = F)
  return(Wk)
}


integralQRdivGdLambda <- function(myZ,myTime,Qvals,dLambda,KMC){
  res <- dplyr::filter(Qvals,myZ==discZ,Ck<=myTime) %>%
    dplyr::mutate(G=KMC(Ck),dLambda=dLambda(Ck)) %>%
    dplyr::summarize(res=sum(QZiCk*dLambda/G)) %>%
    dplyr::select(res) %>%
    unlist(use.names = F)
  return(res)

}



#' @importFrom survival coxph survfit Surv
estimate_competing_risk_model <- function(df){
  model1 <- survival::coxph(survival::Surv(time, status == 1) ~ Z, data = df, method = "breslow")
  b1 <- model1$coefficients["Z"]
  haz1 <- survival::survfit(model1,newdata = data.frame(Z=0))



  model2 <- survival::coxph(Surv(time, status == 2) ~ Z, data = df, method = "breslow")
  b2 <- model2$coefficients["Z"]
  haz2 <- survival::survfit(model2,newdata = data.frame(Z=0))

  return (tibble::tibble(time= haz1$time, b1 = b1, b2 = b2, cumhaz1 = haz1$cumhaz,cumhaz2 = haz2$cumhaz))

}

estimate_cox_model <- function(df){
  return(survival::coxph(survival::Surv(time, status != 0) ~ Z, data = df, method = "breslow"))
}





