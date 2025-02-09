 #install.packages(’survminer’)
 #install.packages(’survivalnma’)
 library(survminer)
 library(survivalnma)


fit <- survfit(Surv(t1_censored[1:(n_AB/2)],delta1[1:(n_AB/2)]) ~ 1)
 d <- data.frame(time = fit$time,
 n.risk = fit$n.risk,
 n.event = fit$n.event
 )


 write.table(d,’Your pathway//s1t1.txt’)
 fit <- survfit(Surv(t1_censored[(n_AB/2+1):n_AB],delta1[(n_AB/2+1):n_AB]) ~ 1)
 d <- data.frame(time = fit$time,
 n.risk = fit$n.risk,
 n.event = fit$n.event
 )
 write.table(d,’Your pathway//s1t2.txt’)
 fit <- survfit(Surv(t2_censored[1:(n_AC/2)],delta2[1:(n_AC/2)]) ~ 1)
 d <- data.frame(time = fit$time,
 n.risk = fit$n.risk,
 n.event = fit$n.event
 )
 write.table(d,’Your pathway//s2t1.txt’)
 fit <- survfit(Surv(t2_censored[(n_AC/2+1):n_AC],delta2[(n_AC/2+1):n_AC]) ~ 1)
 d <- data.frame(time = fit$time,
 n.risk = fit$n.risk, n.event = fit$n.event
 write.table(d,’Your pathway//s2t3.txt’)

  nma_df <- data.frame(
 stringsAsFactors = FALSE,
 "treatment" = c("A", "B","A", "C"),
 "study" = c("Study␣1", "Study␣1", "Study␣2", "Study␣2"),
 "baseline" = c("A", "A", "A", "A"),
 "filepath" = c("your pathway//s1t1.txt",
 "Your pathway//s1t2.txt",
 "Your pathway//s2t1.txt",
 "Your pathway//s2t3.txt"))
 fit_wbl <- survnma(nma_df, model="weibull",bugs.directory =
 "Directory of Winbugs//WinBUGS14",min_time_change = 0.000001)
 #fit_wbl$fit$mean stores the estimates of parameters;
 #fit_wbl$fit$sims.matrix stores the posterior samples
