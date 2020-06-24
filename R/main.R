#' ipcw_estimator: Calculate the ipcw estimator proposed in Li, Goldberg, \& Zheng
#'
#' See details in the Method section
#'
#' @param df  A dataframe with columns obseved time denoted by 'time', covariate 'Z'
#' 'tau1' and 'tau2' which are the lower and upper times denoted by S and t in the paper,
#' and the cometing risk 'status' which is 0 for censoring, 1 for the event of interest, and 2 for the competing risk=status,
#'
#' @return A tibble with the estimate of the inverse-porbablity estimator for functions m_1 to m_6 discussed in the paper.
#' @import tibble
#' @importFrom stats approxfun
#' @importFrom magrittr %>%
#' @export
#'
ipcw_estimator <- function(df){
  KMlist <- KM_Censoring(df)
  KMC <- KMlist$KMC
  n <- nrow(df)
  res <- tibble(type =1:6, ipcw_est= 0)

  for(type in 1:6){
    m1 <- calc_event(df, type = type)

    df <- dplyr::mutate(df,G=KMC(time),m=m1,val=m/G)

    if (type == 3 ) {
      df <- dplyr::mutate(df, t1lesst2= (tau1<=tau2), val =  t1lesst2 + val)
    } else if (type ==6){
      df <- dplyr::mutate(df, t1lesst2= (tau1<=tau2), val =  (1-t1lesst2) + val)
    }
    res$ipcw_est[type] <- sum(df$val)/nrow(df)
  }
  return(res)
}


#' dr_estimator: Calculate the doubly robust estimator proposed in Li, Goldberg, \& Zheng
#'
#' See details in the Method section
#'
#' @param df  A dataframe with columns obseved time denoted by 'time', covariate 'Z'
#' 'tau1' and 'tau2' which are the lower and upper times denoted by S and t in the paper,
#' and the cometing risk 'status' which is 0 for censoring, 1 for the event of interest, and 2 for the competing risk=status,
#' @param discretize_Z  If 'TRUE', discretizes the covariate 'Z'. This can ease the calculation burden
#'  for large dataframe 'df'. Default 'FALSE'. estimate the variance of the estimator. Default 'FALSE'
#'
#' @return A tibble with the estimate of the inverse-porbablity and the doubly roubst estimator for functions m_1 to m_6 discussed in the paper.
#' @import purrr
#' @export
dr_estimator <- function(df, discretize_Z = F){

  est <- ipcw_estimator(df)
  est$dr_est <- 0
  n <- nrow(df)
  # Compute a CIF for every value of Z
  competing_risk_model <- estimate_competing_risk_model(df)
  #cox_model <- estimate_cox_model(df)


  if(discretize_Z){
    df <- dplyr::mutate(df,discZ=discretize_Z_fun(Z))
  } else {
    df <- dplyr::mutate(df,discZ =Z )
  }

  discretize_df <- dplyr::group_by(df,discZ) %>%
    summarise(n=n(),tau1=mean(tau1),tau2=mean(tau2))
  discretize_df <- dplyr::mutate(discretize_df,
                                 CIF1= purrr::map(discZ, .f = CIF, competing_risk_model, event = 1),
                                 CIF2= purrr::map(discZ, .f = CIF, competing_risk_model, event = 2),
                                 S = purrr::map(discZ,.f = SgivenZ, competing_risk_model)
  )



  # Compute the censoring process
  # KMC - kaplan-Meier, R (at risk process), N (counting process)
  # Lambda - cumulative hazard


  KMlist <- KM_Censoring(df)
  Ck <- dplyr::filter(df,status==0) %>%
    dplyr::select(C=time) %>%
    dplyr::arrange(C) %>%
    unlist(use.names = F) %>%
    unique()

  KMC <- KMlist$KMC
  R <- KMlist$R
  N <- KMlist$N
  Lambda <- KMlist$Lambda
  eps <- 1e-10


  Qvals0 <- expand.grid(discretize_df$discZ,Ck) %>%
    tibble::as_tibble() %>%
    dplyr::rename(Z=Var1,Ck=Var2)
  Qvals0 <- dplyr::left_join(Qvals0,discretize_df,by=c("Z"="discZ"))



  Qvals0 <- dplyr::mutate(Qvals0, t1lesst2= (tau1<=tau2),
                          lower12=pmin(pmax(Ck-eps,tau1),tau2),
                          upper= tau2,
                          lower45=pmin(Ck-eps, tau2),
                          lower3.lh = pmin(Ck-eps, tau2),
                          upper45.lh = pmin(tau1, tau2),
                          lower45.lh = pmin(Ck-eps, pmin(tau1, tau2))
  )

  for(type in 1:6){

    switch (as.character(type),
            "1" = {Qvals <- dplyr::mutate(Qvals0,
                                          QZiCk= purrr::pmap_dbl( list(lower12,upper,S,CIF1,Ck,t1lesst2),
                                                                  .f=function(lower,upper,S,CIF1,Ck,t1lesst2){t1lesst2* (CIF1(upper)-CIF1(lower) )/(S(Ck-eps))}))

            },
            "2" = {Qvals <- dplyr::mutate(Qvals0,
                                          QZiCk= purrr::pmap_dbl( list(lower12,upper,S,CIF2,Ck,t1lesst2),
                                                                  .f=function(lower,upper,S,CIF2,Ck,t1lesst2){t1lesst2* (CIF2(upper)-CIF2(lower) )/(S(Ck-eps))}))
            },
            "3" = {Qvals <- dplyr::mutate(Qvals0,
                                          QZiCk= purrr::pmap_dbl( list(lower3.lh,upper,S,Ck,t1lesst2),
                                                                  .f=function(lower3,upper,S,Ck,t1lesst2){-t1lesst2*(S(lower3) - S(upper))/(S(Ck-eps))}))
            },
            "4" = {Qvals <- dplyr::mutate(Qvals0,
                                          QZiCk= purrr::pmap_dbl( list(lower45.lh,upper45.lh,S,CIF1,Ck,t1lesst2),
                                                                  .f=function(lower45,upper,S,CIF1,Ck,t1lesst2){ (CIF1(upper)-CIF1(lower45) )/(S(Ck-eps))}))

            },
            "5" = {Qvals <- dplyr::mutate(Qvals0,
                                          QZiCk= purrr::pmap_dbl( list(lower45.lh,upper45.lh,S,CIF2,Ck,t1lesst2),
                                                                  .f=function(lower45,upper,S,CIF2,Ck,t1lesst2){(CIF2(upper)-CIF2(lower45) )/(S(Ck-eps))}))

            },
            "6" = {Qvals <- dplyr::mutate(Qvals0,
                                          QZiCk= purrr::pmap_dbl( list(lower3.lh,upper,S,Ck,t1lesst2),
                                                                  .f=function(lower3,upper,S,Ck,t1lesst2){-(1-t1lesst2)*(S(lower3) - S(upper))/(S(Ck-eps))})
            )

            }
    )


    Qvals <- dplyr::rename(Qvals,discZ=Z)


    Cvals <- tibble::tibble(Ck=Ck,Rk=R(Ck),dNk=N(Ck)-N(Ck-eps),Gk=KMC(Ck))
    Cvals <- dplyr::mutate(Cvals,
                           sumQZiCk_RiCk=purrr::map_dbl(Ck,.f=sumQZiCk_RiCk,df=df,Qvals=Qvals))


    df0 <- dplyr::filter(df,status==0) %>%
      dplyr::select(C=time,discZ)
    df0 <- dplyr::left_join(df0,Qvals,by=c("C"="Ck","discZ"))
    df0 <- dplyr::left_join(df0,Cvals,by=c("C"="Ck"))
    df0 <- dplyr::mutate(df0,QdivG=QZiCk/Gk,sumQRdivRG=sumQZiCk_RiCk/ (Rk*Gk))


    term1 <- sum(df0$QdivG, na.rm=TRUE)/n

    dLambda <- function(t){Lambda(t)-Lambda(t-1e-8)}

    df1 <- dplyr::mutate(df,intQRdivGdLambda =
                           purrr::map2_dbl(discZ,time,.f = integralQRdivGdLambda,
                                           Qvals=Qvals,dLambda=dLambda,KMC=KMC ))


    term2 <- mean(df1$intQRdivGdLambda, na.rm=TRUE)




    est$dr_est[type] <- min(max(est$ipcw_est[type]+term1-term2,0),1)
  }
  return(est)
}

#' summary_statistics: Calculate summary statistics TPP, FPP, FNP, TNP, Sensitivity, and Specificity
#'
#' See details in the Method section
#'
#' @param est  The result of either the ipcw_estimator or the dr_estimator
#' @return A tibble with the estimate of the inverse-porbablity and the doubly roubst estimator for functions m_1 to m_6 discussed in the paper.
#' @export
summary_statistics <- function(est){
  f <- function(v){
    stat <- numeric(6)
    stat[1] <- v[1] #TPP
    stat[2] <- v[2]+v[3] # FPP
    stat[3] <- v[4] # FNP
    stat[4] <- v[5]+v[6] #TNP
    stat[5] <- stat[1]/(stat[1]+stat[3]) #Sensitivity
    stat[6] <- stat[4]/(stat[2]+stat[4])# Specificity
    return(stat)
  }
  res <- tibble::tibble(statistic = c("TPP","FNP","FNP","TNP","Sensitivity","Specificity")) %>%
    dplyr::mutate(ipcw = f(est$ipcw_est))
  if("dr_est" %in% colnames(est)){
    res <- dplyr::mutate(res,dr= f(est$dr_est))
  }
  return(res)

}


if(getRversion() >= "2.15.1")
  utils::globalVariables(c("m","find.basesurv","b1","cumhaz1","b2","cumhaz2","lam1","S",
                           "lam1prodS","lam2","lam2prodS",
                           'discZ','Z','val','G','n.risk','n.event',"tau1","tau2","status",
                           "time","C","Var1","Var2","lower12","upper","CIF1","t1lesst2",
                           "CIF2","lower3.lh","lower45.lh","QZiCk","Gk","Rk","ipc_est",
                           "Ck", "QZiCk_RiCk", "RiCk", "Wk", "Wk_inner", "dNk","surv", "upper45.lh"))
