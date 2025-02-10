
library(DescTools)
library(expint)
library(survival)
library(flexsurv)
#ns:number of studies
#n_cov: number of covariates
#nt: number of treatments
#n_list: list of sample size.n_list[[i]][[j]] contains the sample size of treatment j in study i.
#t_list: list of t data.For study i with IPD, t_list[[i]][[j]] contains the vector of t of treatment j; For study i with AgD, t_list[[i]][[j]] contains the summary statistics of treatment j in study i
#X_list: list of covariates. For study i with IPD, X_list[[i]][[j]] contains the design matrix of treatment j; For study i with AgD, X_list[[i]][[j]] contains the simulated X of treatment j in study i
#delta_list: list of the truncation information of studies with IPD. delta_list[[i]][[j]] contains the truncation information of treatment j in study i
#n_simu_list: list of number of simulated sample size. Only for study with AgD. n_simu_list[[i]][[j]] contains the simulated sample size for treatment j in study i.
#var_list: for gompertz only; the estimated variance of t for AgD.var_list[[i]][[j]] contains the variance of t of treatment j in study i
#tr_list: list of the treatments in each study; tr_list[[i]] contains the vector of treatments in study i
#st_list: list of the studies containing treatment; st_list[[i]] contains the vector of studies containing treatment i
#sd_pr: sd in the proposal normal distribution
#eta_fre: estimated scale parameter by IPD; used in proposal normal distribution
#n_ipd:number of studies with IPD
#n_agd:number of studies with AgD
#type:weibull_median;weibull_mean;gompertz
#Note: the baseline treatment should be treatment 1;the treatments should be placed in order for tr_list; The study with IPD should be labelled first
mlnmr_survival<-function(ns,n_ipd,n_agd,n_cov,nt,n_list,t_list,X_list,n_simu_list,delta_list,var_list,tr_list,st_list,burnin=20000,nloop=10000,sd_pr=0.1,type,mu_init=rep(8,ns),eta_init=rep(1.2,ns),gamma_init=c(0,rep(-1,nt-1)),beta1_init=rep(0,n_cov),beta2_init=rep(0,n_cov)){
  
  switch(type,
         weibull_median={
           
           f1_ind<-function(n,mu,x,beta1,beta2,gammax,eta,t){
             indx_list<-sapply(1:n,function(i){return((mu*exp(x[i,]%*%(beta1+beta2)+gammax)*t[i]^eta))})
             indx<-sum(indx_list)
             return(indx)
           }
           
           f1_sum<-function(n,n_simu,mu,x_simu,beta1,beta2,gammax,eta,t_sum){
             me_list<-sapply(1:n_simu,function(i){return(mu^(-1/eta)*exp(-(x_simu[i,]%*%(beta1+beta2)/eta+gammax/eta))*(log(2)^(1/eta)))})
             me<-mean(me_list)
             var_simu_list<-sapply(1:n_simu,function(i){return(dweibull(me,eta,1/mu^(1/eta)/exp((x_simu[i,]%*%(beta1+beta2)+gammax)/eta)))})
             f1<-mean(var_simu_list)
             f2<-(t_sum-me)^2*4*n*mean(var_simu_list)^2
             
             return(c(f1,f2))
           }
           
           t<-c()
           delta<-c()
           X<-c()
           tr<-c()
           
           for (i in 1:n_ipd){
             for (j in tr_list[[i]]){
               
               t<-c(t,t_list[[i]][[j]])
               delta<-c(delta,delta_list[[i]][[j]])
               X<-rbind(X,X_list[[i]][[j]])
               tr<-c(tr,rep(j,n_list[[i]][[j]]))
               
               
             }
           }
           
           
           
           fre_res<-survreg(Surv(t,delta)~X*tr,dist='weibull')
           eta_fre<-1/fre_res$scale
           
           
           mu_mc<-mu_init
           eta_mc<-eta_init
           gamma_mc<-gamma_init
           
           beta1_mc<-beta1_init
           beta2_mc<-beta2_init
           
           mu_tab<-matrix(0,nrow=(burnin+nloop),ncol=ns)
           eta_tab<-matrix(0,nrow=(burnin+nloop),ncol=ns)
           gamma_tab<-matrix(0,nrow=(burnin+nloop),ncol=nt)
           beta1_tab<-matrix(0,nrow=(burnin+nloop),ncol=n_cov)
           beta2_tab<-matrix(0,nrow=(burnin+nloop),ncol=n_cov)
           
           for (loop in 1:(burnin+nloop)){
             
             #mu:ipd
             
             mu_mc_tp<-sapply(1:n_ipd,function(i){
               mu_new<-rnorm(1,mu_mc[i],sd_pr)
               
               ipd_v<-0
               if (tr_list[[i]][1]==1){
                 ipd_v<-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])-f1_ind(n_list[[i]][[1]],mu_new,X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                 for (j in  tr_list[[i]][-1]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                 }
               }else{
                 for (j in  tr_list[[i]]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                 }
               }
               
               delta<-c()
               for (j in  tr_list[[i]]){
                 delta<-c(delta,delta_list[[i]][[j]])
               }
               
               r<-exp(ipd_v)*prod((mu_new/mu_mc[i])^(delta))
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,mu_new,mu_mc[i]))
             })
             
             mu_mc[1:n_ipd]<-mu_mc_tp
             
             #mu:agd
             
             mu_mc_tp<-sapply((n_ipd+1):ns,function(i){
               mu_new<-rnorm(1,mu_mc[i],sd_pr)
               
               agd_v_sum<-0
               agd_v_prod<-1
               
               if (tr_list[[i]][1]==1){
                 MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_new,X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                 MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                 agd_v_sum<-MH_mc[2]-MH_new[2]
                 agd_v_prod<-MH_new[1]/MH_mc[1]
                 
                 for (j in tr_list[[i]][-1]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
                 
               }else{
                 for (j in tr_list[[i]]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
               }
               
               r<-exp(agd_v_sum)*(agd_v_prod)#*dnorm(mu2_mc,mu_fre,sd_pr)/dnorm(mu2_new,mu_fre,sd_pr)
               r<-ifelse(is.na(r),0,min(r,1))
               u<-runif(1)
               return(ifelse(u<r,mu_new,mu_mc[i]))
             })
             
             mu_mc[(n_ipd+1):ns]<-mu_mc_tp
             
             # eta:ipd
             
             eta_mc_tp<-sapply(1:n_ipd,function(i){
               
               eta_new<-rnorm(1,eta_fre,0.1)
               ipd_v<-0
               if (tr_list[[i]][1]==1){
                 ipd_v<-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_new,t_list[[i]][[1]])
                 
                 for (j in tr_list[[i]][-1]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]])
                 }
                 
               }else{
                 for (j in tr_list[[i]]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]])
                 }
               }
               
               tc<-c()
               for (j in tr_list[[i]]){
                 tc<-c(tc,t_list[[i]][[j]])
               }
               
               delta<-c()
               for (j in tr_list[[i]]){
                 delta<-c(delta,delta_list[[i]][[j]])
               }
               
               r<-exp(ipd_v)*prod((eta_new/eta_mc[i]*(tc)^(eta_new-eta_mc[i]))^delta)*dnorm(eta_mc[i],eta_fre,0.1)/dnorm(eta_new,eta_fre,0.1)
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,eta_new,eta_mc[i]))
               
               
             })
             
             eta_mc[1:n_ipd]<-eta_mc_tp
             
             
             
             
             #eta:agd
             
             
             eta_mc_tp<-sapply((n_ipd+1):ns,function(i){
               
               eta_new<-rnorm(1,eta_fre,0.1)
               agd_v_sum<-0
               agd_v_prod<-1
               if (tr_list[[i]][1]==1){
                 MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_new,t_list[[i]][[1]])
                 MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                 agd_v_sum<-MH_mc[2]-MH_new[2]
                 agd_v_prod<-MH_new[1]/MH_mc[1]
                 
                 for (j in tr_list[[i]][-1]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
                 
               }else{
                 for (j in tr_list[[i]]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
               }
               
               r<-exp(agd_v_sum)*(agd_v_prod)*(gamma(1+2/eta_mc[i])-gamma(1+1/eta_mc[i])^2)/(gamma(1+2/eta_new)-gamma(1+1/eta_new)^2)*dnorm(eta_mc[i],eta_fre,0.1)/dnorm(eta_new,eta_fre,0.1)
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,eta_new,eta_mc[i]))
               
               
             })
             
             eta_mc[(n_ipd+1):ns]<-eta_mc_tp
             
             #gamma
             
             gamma_mc_tp<-sapply(2:nt,function(j){
               gamma_new<-rnorm(1,gamma_mc[j],0.01)
               
               ipd_s<-st_list[[j]][which(st_list[[j]]<=n_ipd)]
               agd_s<-st_list[[j]][which(st_list[[j]]>n_ipd)]
               
               ipd_r<-1
               if (length(ipd_s)>0){
                 for (i in ipd_s){
                   j<-j
                   ipd_v<-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_new,eta_mc[i],t_list[[i]][[j]])
                   ipd_r<-ipd_r*exp(ipd_v)*prod(exp(gamma_new-gamma_mc[j])^delta_list[[i]][[j]])
                 }
               }
               
               agd_r<-1
               if (length(agd_s)>0){
                 for (i in agd_s){
                   j<-j
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_new,eta_mc[i],t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   agd_v_sum<-MH_mc[2]-MH_new[2]
                   agd_v_prod<-MH_new[1]/MH_mc[1]
                   agd_r<-agd_r*exp(agd_v_sum)*agd_v_prod
                 }
               }
               
               r<-agd_r*ipd_r
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,gamma_new,gamma_mc[j]))
             })
             
             gamma_mc[2:(nt)]<-gamma_mc_tp
             
             
             #beta1_mc
             
             for (l in 1:length(beta1_mc)){
               beta1_new<-beta1_mc
               beta1_new[l]<-rnorm(1,beta1_mc[l],0.1*sd_pr)
               
               
               beta1_ipd_tp<-sapply(1:n_ipd,function(i){
                 ipd_v<-0
                 if (tr_list[[i]][1]==1){
                   
                   ipd_v<-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_new,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                   for (j in tr_list[[i]][-1]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                 }
                 
                 X<-c()
                 for (j in tr_list[[i]]){
                   X<-rbind(X,X_list[[i]][[j]])
                 }
                 
                 delta<-c()
                 for (j in tr_list[[i]]){
                   delta<-c(delta,delta_list[[i]][[j]])
                 }
                 
                 return(exp(ipd_v)*prod((exp(X%*%(beta1_new-beta1_mc)))^(delta)))
               })
               
               beta1_agd_tp<-sapply((n_ipd+1):ns,function(i){
                 
                 agd_v_sum<-0
                 agd_v_prod<-1
                 if (tr_list[[i]][1]==1){
                   
                   MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_new,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                   MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                   
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                   for (j in tr_list[[i]][-1]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }
                 
                 
                 return(c(agd_v_sum,agd_v_prod))
               })
               
               
               agd_v_sum<-sum(beta1_agd_tp[1,])
               agd_v_prod<-prod(beta1_agd_tp[2,])
               
               
               r<-prod(beta1_ipd_tp)*exp(agd_v_sum)*agd_v_prod
               r<-ifelse(is.na(r),0,min(r,1))
               beta1_mc[l]<-ifelse(runif(1)<r,beta1_new[l],beta1_mc[l])
               
             }
             
             #beta2_mc
             
             for (l in 1:length(beta2_mc)){
               beta2_new<-beta2_mc
               beta2_new[l]<-rnorm(1,beta2_mc[l],0.1*sd_pr)
               
               
               beta2_ipd_tp<-sapply(1:n_ipd,function(i){
                 ipd_v<-0
                 if (tr_list[[i]][1]==1){
                   
                   for (j in tr_list[[i]][-1]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                   
                   delta<-c()
                   for (j in tr_list[[i]][-1]){
                     delta<-c(delta,delta_list[[i]][[j]])
                   }
                   
                   X<-c()
                   for (j in tr_list[[i]][-1]){
                     X<-rbind(X,X_list[[i]][[j]])
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                   
                   delta<-c()
                   for (j in tr_list[[i]]){
                     delta<-c(delta,delta_list[[i]][[j]])
                   }
                   
                   X<-c()
                   for (j in tr_list[[i]]){
                     X<-rbind(X,X_list[[i]][[j]])
                   }
                 }
                 
                 
                 
                 
                 
                 return(exp(ipd_v)*prod((exp(X%*%(beta2_new-beta2_mc)))^(delta)))
               })
               
               beta2_agd_tp<-sapply((n_ipd+1):ns,function(i){
                 
                 agd_v_sum<-0
                 agd_v_prod<-1
                 if (tr_list[[i]][1]==1){
                   
                   for (j in tr_list[[i]][-1]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }
                 
                 
                 return(c(agd_v_sum,agd_v_prod))
               })
               
               
               agd_v_sum<-sum(beta2_agd_tp[1,])
               agd_v_prod<-prod(beta2_agd_tp[2,])
               
               
               r<-prod(beta2_ipd_tp)*exp(agd_v_sum)*agd_v_prod
               r<-ifelse(is.na(r),0,min(r,1))
               beta2_mc[l]<-ifelse(runif(1)<r,beta2_new[l],beta2_mc[l])
               
             }
             
             mu_tab[loop,]<-mu_mc
             eta_tab[loop,]<-eta_mc
             gamma_tab[loop,]<-gamma_mc
             beta1_tab[loop,]<-beta1_mc
             beta2_tab[loop,]<-beta2_mc
             
           }
           
           
           mu_est<-apply(mu_tab[(burnin+1):(burnin+nloop),],2,mean)
           eta_est<-apply(eta_tab[(burnin+1):(burnin+nloop),],2,mean)
           gamma_est<-apply(gamma_tab[(burnin+1):(burnin+nloop),],2,mean)
           beta1_est<-apply(beta1_tab[(burnin+1):(burnin+nloop),],2,mean)
           beta2_est<-apply(beta2_tab[(burnin+1):(burnin+nloop),],2,mean)
           
           
           
           return(list(mu_est,eta_est,gamma_est,beta1_est,beta2_est,mu_tab[(burnin+1):(burnin+nloop),],eta_tab[(burnin+1):(burnin+nloop),],gamma_tab[(burnin+1):(burnin+nloop),],beta1_tab[(burnin+1):(burnin+nloop),],beta2_tab[(burnin+1):(burnin+nloop),]))
         },weibull_mean={
           
           f1_ind<-function(n,mu,x,beta1,beta2,gammax,eta,t){
             indx_list<-sapply(1:n,function(i){return((mu*exp(x[i,]%*%(beta1+beta2)+gammax)*t[i]^eta))})
             indx<-sum(indx_list)
             return(indx)
           }
           
           f1_sum<-function(n,n_simu,mu,x_simu,beta1,beta2,gammax,eta,t_sum){
             var_simu_list<-sapply(1:n,function(i){return(exp(-2*(x_simu[i,]%*%(beta1+beta2)/eta+gammax/eta)))})
             mean_simu_list<-sapply(1:n,function(i){return(exp(-(x_simu[i,]%*%(beta1+beta2)/eta+gammax/eta)))})
             
             f1<-mu^(2/eta)/sqrt(mean(var_simu_list))
             f2<-((t_sum-mean(mean_simu_list)/(mu^(1/eta))*gamma(1+1/eta))^2)*mu^(2/eta)/(2/n*mean(var_simu_list)*(gamma(1+2/eta)-gamma(1+1/eta)^2))
             
             return(c(f1,f2))
           } 
           
           t<-c()
           delta<-c()
           X<-c()
           tr<-c()
           
           for (i in 1:n_ipd){
             for (j in tr_list[[i]]){
               
               t<-c(t,t_list[[i]][[j]])
               delta<-c(delta,delta_list[[i]][[j]])
               X<-rbind(X,X_list[[i]][[j]])
               tr<-c(tr,rep(j,n_list[[i]][[j]]))
               
               
             }
           }
           
           
           
           fre_res<-survreg(Surv(t,delta)~X*tr,dist='weibull')
           eta_fre<-1/fre_res$scale
           
           mu_mc<-mu_init
           eta_mc<-eta_init
           gamma_mc<-gamma_init
           
           beta1_mc<-beta1_init
           beta2_mc<-beta2_init
           
           mu_tab<-matrix(0,nrow=(burnin+nloop),ncol=ns)
           eta_tab<-matrix(0,nrow=(burnin+nloop),ncol=ns)
           gamma_tab<-matrix(0,nrow=(burnin+nloop),ncol=nt)
           beta1_tab<-matrix(0,nrow=(burnin+nloop),ncol=n_cov)
           beta2_tab<-matrix(0,nrow=(burnin+nloop),ncol=n_cov)
           
           for (loop in 1:(burnin+nloop)){
             
             #mu:ipd
             
             mu_mc_tp<-sapply(1:n_ipd,function(i){
               mu_new<-rnorm(1,mu_mc[i],sd_pr)
               
               ipd_v<-0
               if (tr_list[[i]][1]==1){
                 ipd_v<-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])-f1_ind(n_list[[i]][[1]],mu_new,X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                 for (j in  tr_list[[i]][-1]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                 }
               }else{
                 for (j in  tr_list[[i]]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                 }
               }
               
               delta<-c()
               for (j in  tr_list[[i]]){
                 delta<-c(delta,delta_list[[i]][[j]])
               }
               
               r<-exp(ipd_v)*prod((mu_new/mu_mc[i])^(delta))
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,mu_new,mu_mc[i]))
             })
             
             mu_mc[1:n_ipd]<-mu_mc_tp
             
             #mu:agd
             
             mu_mc_tp<-sapply((n_ipd+1):ns,function(i){
               mu_new<-rnorm(1,mu_mc[i],sd_pr)
               
               agd_v_sum<-0
               agd_v_prod<-1
               
               if (tr_list[[i]][1]==1){
                 MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_new,X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                 MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                 agd_v_sum<-MH_mc[2]-MH_new[2]
                 agd_v_prod<-MH_new[1]/MH_mc[1]
                 
                 for (j in tr_list[[i]][-1]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
                 
               }else{
                 for (j in tr_list[[i]]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
               }
               
               r<-exp(agd_v_sum)*(agd_v_prod)#*dnorm(mu2_mc,mu_fre,sd_pr)/dnorm(mu2_new,mu_fre,sd_pr)
               r<-ifelse(is.na(r),0,min(r,1))
               u<-runif(1)
               return(ifelse(u<r,mu_new,mu_mc[i]))
             })
             
             mu_mc[(n_ipd+1):ns]<-mu_mc_tp
             
             # eta:ipd
             
             eta_mc_tp<-sapply(1:n_ipd,function(i){
               
               eta_new<-rnorm(1,eta_fre,0.1)
               ipd_v<-0
               if (tr_list[[i]][1]==1){
                 ipd_v<-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_new,t_list[[i]][[1]])
                 
                 for (j in tr_list[[i]][-1]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]])
                 }
                 
               }else{
                 for (j in tr_list[[i]]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]])
                 }
               }
               
               tc<-c()
               for (j in tr_list[[i]]){
                 tc<-c(tc,t_list[[i]][[j]])
               }
               
               delta<-c()
               for (j in tr_list[[i]]){
                 delta<-c(delta,delta_list[[i]][[j]])
               }
               
               r<-exp(ipd_v)*prod((eta_new/eta_mc[i]*(tc)^(eta_new-eta_mc[i]))^delta)*dnorm(eta_mc[i],eta_fre,0.1)/dnorm(eta_new,eta_fre,0.1)
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,eta_new,eta_mc[i]))
               
               
             })
             
             eta_mc[1:n_ipd]<-eta_mc_tp
             
             
             
             
             #eta:agd
             
             
             eta_mc_tp<-sapply((n_ipd+1):ns,function(i){
               
               eta_new<-rnorm(1,eta_fre,0.1)
               agd_v_sum<-0
               agd_v_prod<-1
               if (tr_list[[i]][1]==1){
                 MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_new,t_list[[i]][[1]])
                 MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                 agd_v_sum<-MH_mc[2]-MH_new[2]
                 agd_v_prod<-MH_new[1]/MH_mc[1]
                 
                 for (j in tr_list[[i]][-1]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
                 
               }else{
                 for (j in tr_list[[i]]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
               }
               
               r<-exp(agd_v_sum)*(agd_v_prod)*(gamma(1+2/eta_mc[i])-gamma(1+1/eta_mc[i])^2)/(gamma(1+2/eta_new)-gamma(1+1/eta_new)^2)*dnorm(eta_mc[i],eta_fre,0.1)/dnorm(eta_new,eta_fre,0.1)
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,eta_new,eta_mc[i]))
               
               
             })
             
             eta_mc[(n_ipd+1):ns]<-eta_mc_tp
             
             #gamma
             
             gamma_mc_tp<-sapply(2:nt,function(j){
               gamma_new<-rnorm(1,gamma_mc[j],0.01)
               
               ipd_s<-st_list[[j]][which(st_list[[j]]<=n_ipd)]
               agd_s<-st_list[[j]][which(st_list[[j]]>n_ipd)]
               
               ipd_r<-1
               if (length(ipd_s)>0){
                 for (i in ipd_s){
                   j<-j
                   ipd_v<-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_new,eta_mc[i],t_list[[i]][[j]])
                   ipd_r<-ipd_r*exp(ipd_v)*prod(exp(gamma_new-gamma_mc[j])^delta_list[[i]][[j]])
                 }
               }
               
               agd_r<-1
               if (length(agd_s)>0){
                 for (i in agd_s){
                   j<-j
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_new,eta_mc[i],t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   agd_v_sum<-MH_mc[2]-MH_new[2]
                   agd_v_prod<-MH_new[1]/MH_mc[1]
                   agd_r<-agd_r*exp(agd_v_sum)*agd_v_prod
                 }
               }
               
               r<-agd_r*ipd_r
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,gamma_new,gamma_mc[j]))
             })
             
             gamma_mc[2:(nt)]<-gamma_mc_tp
             
             
             #beta1_mc
             
             for (l in 1:length(beta1_mc)){
               beta1_new<-beta1_mc
               beta1_new[l]<-rnorm(1,beta1_mc[l],0.1*sd_pr)
               
               
               beta1_ipd_tp<-sapply(1:n_ipd,function(i){
                 ipd_v<-0
                 if (tr_list[[i]][1]==1){
                   
                   ipd_v<-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_new,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                   for (j in tr_list[[i]][-1]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                 }
                 
                 X<-c()
                 for (j in tr_list[[i]]){
                   X<-rbind(X,X_list[[i]][[j]])
                 }
                 
                 delta<-c()
                 for (j in tr_list[[i]]){
                   delta<-c(delta,delta_list[[i]][[j]])
                 }
                 
                 return(exp(ipd_v)*prod((exp(X%*%(beta1_new-beta1_mc)))^(delta)))
               })
               
               beta1_agd_tp<-sapply((n_ipd+1):ns,function(i){
                 
                 agd_v_sum<-0
                 agd_v_prod<-1
                 if (tr_list[[i]][1]==1){
                   
                   MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_new,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                   MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                   
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                   for (j in tr_list[[i]][-1]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }
                 
                 
                 return(c(agd_v_sum,agd_v_prod))
               })
               
               
               agd_v_sum<-sum(beta1_agd_tp[1,])
               agd_v_prod<-prod(beta1_agd_tp[2,])
               
               
               r<-prod(beta1_ipd_tp)*exp(agd_v_sum)*agd_v_prod
               r<-ifelse(is.na(r),0,min(r,1))
               beta1_mc[l]<-ifelse(runif(1)<r,beta1_new[l],beta1_mc[l])
               
             }
             
             #beta2_mc
             
             for (l in 1:length(beta2_mc)){
               beta2_new<-beta2_mc
               beta2_new[l]<-rnorm(1,beta2_mc[l],0.1*sd_pr)
               
               
               beta2_ipd_tp<-sapply(1:n_ipd,function(i){
                 ipd_v<-0
                 if (tr_list[[i]][1]==1){
                   
                   for (j in tr_list[[i]][-1]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                   
                   delta<-c()
                   for (j in tr_list[[i]][-1]){
                     delta<-c(delta,delta_list[[i]][[j]])
                   }
                   
                   X<-c()
                   for (j in tr_list[[i]][-1]){
                     X<-rbind(X,X_list[[i]][[j]])
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                   
                   delta<-c()
                   for (j in tr_list[[i]]){
                     delta<-c(delta,delta_list[[i]][[j]])
                   }
                   
                   X<-c()
                   for (j in tr_list[[i]]){
                     X<-rbind(X,X_list[[i]][[j]])
                   }
                 }
                 
                 
                 
                 
                 
                 return(exp(ipd_v)*prod((exp(X%*%(beta2_new-beta2_mc)))^(delta)))
               })
               
               beta2_agd_tp<-sapply((n_ipd+1):ns,function(i){
                 
                 agd_v_sum<-0
                 agd_v_prod<-1
                 if (tr_list[[i]][1]==1){
                   
                   for (j in tr_list[[i]][-1]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }
                 
                 
                 return(c(agd_v_sum,agd_v_prod))
               })
               
               
               agd_v_sum<-sum(beta2_agd_tp[1,])
               agd_v_prod<-prod(beta2_agd_tp[2,])
               
               
               r<-prod(beta2_ipd_tp)*exp(agd_v_sum)*agd_v_prod
               r<-ifelse(is.na(r),0,min(r,1))
               beta2_mc[l]<-ifelse(runif(1)<r,beta2_new[l],beta2_mc[l])
               
             }
             
             mu_tab[loop,]<-mu_mc
             eta_tab[loop,]<-eta_mc
             gamma_tab[loop,]<-gamma_mc
             beta1_tab[loop,]<-beta1_mc
             beta2_tab[loop,]<-beta2_mc
             
           }
           
           
           mu_est<-apply(mu_tab[(burnin+1):(burnin+nloop),],2,mean)
           eta_est<-apply(eta_tab[(burnin+1):(burnin+nloop),],2,mean)
           gamma_est<-apply(gamma_tab[(burnin+1):(burnin+nloop),],2,mean)
           beta1_est<-apply(beta1_tab[(burnin+1):(burnin+nloop),],2,mean)
           beta2_est<-apply(beta2_tab[(burnin+1):(burnin+nloop),],2,mean)
           
           
           
           return(list(mu_est,eta_est,gamma_est,beta1_est,beta2_est,mu_tab[(burnin+1):(burnin+nloop),],eta_tab[(burnin+1):(burnin+nloop),],gamma_tab[(burnin+1):(burnin+nloop),],beta1_tab[(burnin+1):(burnin+nloop),],beta2_tab[(burnin+1):(burnin+nloop),]))
         },gompertz={
           
          
           f1_ind<-function(n,mu,x,beta1,beta2,gammax,eta,t){
             indx_list<-sapply(1:n,function(i){return(mu*exp(x[i,]%*%(beta1+beta2)+gammax)/eta*(exp(eta*t[i])-1))})
             indx<-sum(indx_list)
             return(indx)
           }
           
           f1_sum<-function(n,n_simu,mu,x_simu,beta1,beta2,gammax,eta,t_sum,t_var){
             simu_list<-sapply(1:n_simu,function(i){
               
               shape<-eta;rate<-mu/eta*exp(x_simu[i,]%*%(beta1+beta2)+gammax)
               
               simu_mean<-1/shape*exp(rate)*(-expint_E1(-rate))
               
               return(simu_mean)
             })
             
             
             f1<-1/sqrt(t_var)
             f2<-((t_sum-mean(simu_list))^2)*n/2/t_var
             
             return(c(f1,f2))
           }     
           
           
           t<-c()
           delta<-c()
           X<-c()
           tr<-c()
           
           for (i in 1:n_ipd){
             for (j in tr_list[[i]]){
               
               t<-c(t,t_list[[i]][[j]])
               delta<-c(delta,delta_list[[i]][[j]])
               X<-rbind(X,X_list[[i]][[j]])
               tr<-c(tr,rep(j,n_list[[i]][[j]]))
               
               
             }
           }
           
           fre_res<-flexsurvreg(formula = Surv(t,delta)~tr*X, dist="gompertz")
           eta_fre<-fre_res$coefficients[1]
           
           mu_mc<-mu_init
           eta_mc<-eta_init
           gamma_mc<-gamma_init
           
           beta1_mc<-beta1_init
           beta2_mc<-beta2_init
           
           mu_tab<-matrix(0,nrow=(burnin+nloop),ncol=ns)
           eta_tab<-matrix(0,nrow=(burnin+nloop),ncol=ns)
           gamma_tab<-matrix(0,nrow=(burnin+nloop),ncol=nt)
           beta1_tab<-matrix(0,nrow=(burnin+nloop),ncol=n_cov)
           beta2_tab<-matrix(0,nrow=(burnin+nloop),ncol=n_cov)
           
           for (loop in 1:(burnin+nloop)){
             
             #mu:ipd
             
             mu_mc_tp<-sapply(1:n_ipd,function(i){
               mu_new<-rnorm(1,mu_mc[i],sd_pr)
               
               ipd_v<-0
               if (tr_list[[i]][1]==1){
                 ipd_v<-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])-f1_ind(n_list[[i]][[1]],mu_new,X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                 for (j in tr_list[[i]][-1]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                 }
               }else{
                 for (j in tr_list[[i]]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                 }
               }
               
               delta<-c()
               for (j in tr_list[[i]]){
                 delta<-c(delta,delta_list[[i]][[j]])
               }
               
               r<-exp(ipd_v)*prod((mu_new/mu_mc[i])^(delta))
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,mu_new,mu_mc[i]))
             })
             
             mu_mc[1:n_ipd]<-mu_mc_tp
             
             #mu:agd
             
             mu_mc_tp<-sapply((n_ipd+1):ns,function(i){
               mu_new<-rnorm(1,mu_mc[i],sd_pr)
               
               agd_v_sum<-0
               agd_v_prod<-1
               
               if (tr_list[[i]][1]==1){
                 MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_new,X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]],var_list[[i]][[1]])
                 MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]],var_list[[i]][[1]])
                 agd_v_sum<-MH_mc[2]-MH_new[2]
                 agd_v_prod<-MH_new[1]/MH_mc[1]
                 
                 for (j in tr_list[[i]][-1]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
                 
               }else{
                 for (j in tr_list[[i]]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
               }
               
               r<-exp(agd_v_sum)*(agd_v_prod)#*dnorm(mu2_mc,mu_fre,sd_pr)/dnorm(mu2_new,mu_fre,sd_pr)
               r<-ifelse(is.na(r),0,min(r,1))
               u<-runif(1)
               return(ifelse(u<r,mu_new,mu_mc[i]))
             })
             
             mu_mc[(n_ipd+1):ns]<-mu_mc_tp
             
             # eta:ipd
             
             eta_mc_tp<-sapply(1:n_ipd,function(i){
               
               eta_new<-rnorm(1,eta_fre,0.1)
               ipd_v<-0
               if (tr_list[[i]][1]==1){
                 ipd_v<-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_new,t_list[[i]][[1]])
                 
                 for (j in tr_list[[i]][-1]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]])
                 }
                 
               }else{
                 for (j in tr_list[[i]]){
                   ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]])
                 }
               }
               
               tc<-c()
               for (j in tr_list[[i]]){
                 tc<-c(tc,t_list[[i]][[j]])
               }
               
               delta<-c()
               for (j in tr_list[[i]]){
                 delta<-c(delta,delta_list[[i]][[j]])
               }
               
               r<-exp(ipd_v)*exp(sum((eta_new-eta_mc[i])*(tc)*delta))*dnorm(eta_mc[i],eta_fre,0.1)/dnorm(eta_new,eta_fre,0.1)
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,eta_new,eta_mc[i]))
               
               
             })
             
             eta_mc[1:n_ipd]<-eta_mc_tp
             
             
             
             
             #eta:agd
             
             
             eta_mc_tp<-sapply((n_ipd+1):ns,function(i){
               
               eta_new<-rnorm(1,eta_fre,0.1)
               agd_v_sum<-0
               agd_v_prod<-1
               if (tr_list[[i]][1]==1){
                 MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_new,t_list[[i]][[1]],var_list[[i]][[1]])
                 MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]],var_list[[i]][[1]])
                 agd_v_sum<-MH_mc[2]-MH_new[2]
                 agd_v_prod<-MH_new[1]/MH_mc[1]
                 
                 for (j in tr_list[[i]][-1]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]],var_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
                 
               }else{
                 for (j in tr_list[[i]]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_new,t_list[[i]][[j]],var_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                 }
               }
               
               r<-exp(agd_v_sum)*(agd_v_prod)
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,eta_new,eta_mc[i]))
               
               
             })
             
             eta_mc[(n_ipd+1):ns]<-eta_mc_tp
             
             #gamma
             
             gamma_mc_tp<-sapply(2:nt,function(j){
               gamma_new<-rnorm(1,gamma_mc[j],0.01)
               
               ipd_s<-st_list[[j]][which(st_list[[j]]<=n_ipd)]
               agd_s<-st_list[[j]][which(st_list[[j]]>n_ipd)]
               
               ipd_r<-1
               if (length(ipd_s)>0){
                 for (i in ipd_s){
                   
                   ipd_v<-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_new,eta_mc[i],t_list[[i]][[j]])
                   ipd_r<-ipd_r*exp(ipd_v)*prod(exp(gamma_new-gamma_mc[j])^delta_list[[i]][[j]])
                 }
               }
               
               agd_r<-1
               if (length(agd_s)>0){
                 for (i in agd_s){
                   
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_new,eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                   agd_v_sum<-MH_mc[2]-MH_new[2]
                   agd_v_prod<-MH_new[1]/MH_mc[1]
                   agd_r<-agd_r*exp(agd_v_sum)*agd_v_prod
                 }
               }
               
               r<-agd_r*ipd_r
               return(ifelse(runif(1)<r,gamma_new,gamma_mc[j]))
             })
             
             gamma_mc[2:(nt)]<-gamma_mc_tp
             
             
             #beta1_mc
             
             for (l in 1:length(beta1_mc)){
               beta1_new<-beta1_mc
               beta1_new[l]<-rnorm(1,beta1_mc[l],0.1*sd_pr)
               
               
               beta1_ipd_tp<-sapply(1:n_ipd,function(i){
                 ipd_v<-0
                 if (tr_list[[i]][1]==1){
                   
                   ipd_v<-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])-f1_ind(n_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_new,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]])
                   for (j in tr_list[[i]][-1]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                 }
                 
                 X<-c()
                 for (j in tr_list[[i]]){
                   X<-rbind(X,X_list[[i]][[j]])
                 }
                 
                 delta<-c()
                 for (j in tr_list[[i]]){
                   delta<-c(delta,delta_list[[i]][[j]])
                 }
                 
                 return(exp(ipd_v)*prod((exp(X%*%(beta1_new-beta1_mc)))^(delta)))
               })
               
               beta1_agd_tp<-sapply((n_ipd+1):ns,function(i){
                 
                 agd_v_sum<-0
                 agd_v_prod<-1
                 if (tr_list[[i]][1]==1){
                   
                   MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_new,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]],var_list[[i]][[1]])
                   MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,eta_mc[i],t_list[[i]][[1]],var_list[[i]][[1]])
                   
                   agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                   agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                   for (j in tr_list[[i]][-1]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }
                 
                 
                 return(c(agd_v_sum,agd_v_prod))
               })
               
               
               agd_v_sum<-sum(beta1_agd_tp[1,])
               agd_v_prod<-prod(beta1_agd_tp[2,])
               
               
               r<-prod(beta1_ipd_tp)*exp(agd_v_sum)*agd_v_prod
               r<-ifelse(is.na(r),0,min(r,1))
               beta1_mc[l]<-ifelse(runif(1)<r,beta1_new[l],beta1_mc[l])
               
             }
             
             #beta2_mc
             
             for (l in 1:length(beta2_mc)){
               beta2_new<-beta2_mc
               beta2_new[l]<-rnorm(1,beta2_mc[l],0.1*sd_pr)
               
               
               beta2_ipd_tp<-sapply(1:n_ipd,function(i){
                 ipd_v<-0
                 if (tr_list[[i]][1]==1){
                   
                   for (j in tr_list[[i]][-1]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                   
                   delta<-c()
                   for (j in tr_list[[i]][-1]){
                     delta<-c(delta,delta_list[[i]][[j]])
                   }
                   
                   X<-c()
                   for (j in tr_list[[i]][-1]){
                     X<-rbind(X,X_list[[i]][[j]])
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ipd_v<-ipd_v+f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])-f1_ind(n_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]])
                   }
                   
                   delta<-c()
                   for (j in tr_list[[i]]){
                     delta<-c(delta,delta_list[[i]][[j]])
                   }
                   
                   X<-c()
                   for (j in tr_list[[i]]){
                     X<-rbind(X,X_list[[i]][[j]])
                   }
                 }
                 
                 
                 
                 
                 
                 return(exp(ipd_v)*prod((exp(X%*%(beta2_new-beta2_mc)))^(delta)))
               })
               
               beta2_agd_tp<-sapply((n_ipd+1):ns,function(i){
                 
                 agd_v_sum<-0
                 agd_v_prod<-1
                 if (tr_list[[i]][1]==1){
                   
                   for (j in tr_list[[i]][-1]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_new,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],eta_mc[i],t_list[[i]][[j]],var_list[[i]][[j]])
                     
                     agd_v_sum<-agd_v_sum+MH_mc[2]-MH_new[2]
                     agd_v_prod<-agd_v_prod*MH_new[1]/MH_mc[1]
                     
                   }
                 }
                 
                 
                 return(c(agd_v_sum,agd_v_prod))
               })
               
               
               agd_v_sum<-sum(beta2_agd_tp[1,])
               agd_v_prod<-prod(beta2_agd_tp[2,])
               
               
               r<-prod(beta2_ipd_tp)*exp(agd_v_sum)*agd_v_prod
               r<-ifelse(is.na(r),0,min(r,1))
               beta2_mc[l]<-ifelse(runif(1)<r,beta2_new[l],beta2_mc[l])
               
             }
             
             mu_tab[loop,]<-mu_mc
             eta_tab[loop,]<-eta_mc
             gamma_tab[loop,]<-gamma_mc
             beta1_tab[loop,]<-beta1_mc
             beta2_tab[loop,]<-beta2_mc
             
           }
           
           
           mu_est<-apply(mu_tab[(burnin+1):(burnin+nloop),],2,mean)
           eta_est<-apply(eta_tab[(burnin+1):(burnin+nloop),],2,mean)
           gamma_est<-apply(gamma_tab[(burnin+1):(burnin+nloop),],2,mean)
           beta1_est<-apply(beta1_tab[(burnin+1):(burnin+nloop),],2,mean)
           beta2_est<-apply(beta2_tab[(burnin+1):(burnin+nloop),],2,mean)
           
           
           
           return(list(mu_est,eta_est,gamma_est,beta1_est,beta2_est,mu_tab[(burnin+1):(burnin+nloop),],eta_tab[(burnin+1):(burnin+nloop),],gamma_tab[(burnin+1):(burnin+nloop),],beta1_tab[(burnin+1):(burnin+nloop),],beta2_tab[(burnin+1):(burnin+nloop),]))
           
         }
         
  )
  
  
  
  
}
