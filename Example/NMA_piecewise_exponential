library(MASS)
library(survival)

library(dplyr)
library(gemtc)        
library(gemtcPlus)     
library(ggmcmc) 

fit1 <- survfit(Surv(t1_censored[1:(n_AB/2)],delta1[1:(n_AB/2)]) ~ 1)
fit2 <- survfit(Surv(t1_censored[(n_AB/2+1):(n_AB)],delta1[(n_AB/2+1):(n_AB)]) ~ 1)
fit3 <- survfit(Surv(t2_censored[1:(n_AC/2)],delta2[1:(n_AC/2)]) ~ 1)
fit4 <- survfit(Surv(t2_censored[(n_AC/2+1):(n_AC)],delta2[(n_AC/2+1):(n_AC)]) ~ 1)

study=c(rep('STUDY1',length(fit1$n.event)),rep('STUDY1',length(fit2$n.risk)),rep('STUDY2',length(fit3$n.event)),rep('STUDY2',length(fit4$n.censor)))
treatment=c(rep("A",length(fit1$n.event)),rep('B',length(fit2$n.risk)),rep('A',length(fit3$n.event)),rep('C',length(fit4$n.censor)))
t.start=c(c(0,fit1$time[-length(fit1$time)]),c(0,fit2$time[-length(fit2$time)]),c(0,fit3$time[-length(fit3$time)]),c(0,fit4$time[-length(fit4$time)]))
t.end=c(fit1$time,fit2$time,fit3$time,fit4$time)
n.event=c(fit1$n.event,fit2$n.event,fit3$n.event,fit4$n.event)
n.censored=c(fit1$n.censor,fit2$n.censor,fit3$n.censor,fit4$n.censor)
n.risk=c(fit1$n.risk,fit2$n.risk,fit3$n.risk,fit4$n.risk)
  
  
d <- data.frame(study = relevel(factor(study), ref="STUDY1"),
                  treatment = relevel(factor(treatment), ref="A"),
                  t.start = t.start,
                  t.end = t.end,
                  n.event = n.event,
                  n.censored = n.censored,
                  n.risk = n.risk
)
  
model_plan <- plan_pwe(model.pars = list(cut.pts =  c(0.8)),
                         bth.model = "FE", ref.std = "STUDY1", nma.ref.trt = "A",
                         n.chains = 2,
                         n.iter = 6000,
                         n.burnin = 1000,
                         n.thin = 1)
model_input <- nma_pre_proc(d, model_plan)
model  <- nma_fit(model_input = model_input)
