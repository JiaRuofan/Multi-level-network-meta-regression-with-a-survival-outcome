#ns:number of studies
#n_cov: number of covariates
#nt: number of treatments
#n_list: list of sample size.For IPD, n_list[[i]][[j]][[1]] contains the noncensored sample size of treatment j in study i;n_list[[i]][[j]][[2]] contains the censored sample size of treatment j in study i;For AgD, n_list[[i]][[j]] contains the sample size of treatment j in study i
#t_list: list of transformed t data (log t).For IPD, t_list[[i]][[j]][[1]] contains the vector of noncensored t of treatment j in study i;t_list[[i]][[j]][[2]] contains the vector of censored t of treatment j in study i;For AgD, t_list[[i]][[j]] contains the summary statistics of t of treatment j in study i
#X_list: list of covariates. For IPD, X_list[[i]][[j]][[1]] contains the  design matrix of noncensored sample of treatment j in study i;X_list[[i]][[j]][[2]] contains the design matrix of censored sample of treatment j in study i;For AgD, X_list[[i]][[j]] contains the simulated X of treatment j in study i
#n_simu_list: list of number of simulated sample size. Only for study with AgD. n_simu_list[[i]][[j]] contains the simulated sample size for treatment j in study i.
#tr_list: list of the treatments in each study; tr_list[[i]] contains the vector of treatments in study i
#st_list: list of the studies containing treatment; st_list[[i]] contains the vector of studies containing treatment i
#sd_pr: sd in the proposal normal distribution
#eta_fre: estimated scale parameter by IPD; used in proposal normal distribution
#n_ipd:number of studies with IPD
#n_agd:number of studies with AgD
#type:logn;loglog
#mu_init,eta_init: initial value of mu and eta, a vector of length ns
#gamma_init: initial value of gamma, a vector of length nt with first value being 0. 
#beta1_init,beta2_init: initial value of beta1 and beta2, vectors of length n_cov
#Note: the baseline treatment should be treatment 1, the treatments should be placed in order for tr_list
mlnmr_survival_trans<-function(ns,n_ipd,n_agd,n_cov,nt,n_list,t_list,X_list,n_simu_list,tr_list,st_list,burnin=20000,nloop=10000,sd_pr=0.1,type,mu_init=rep(log(8),ns),sigma_init=rep(1.2,ns),gamma_init=c(0,rep(-1,nt-1)),beta1_init=rep(0,n_cov),beta2_init=rep(0,n_cov)){
  
  switch(type,
         logn={
           
           
           f1_ind_ncen<-function(n,mu,x,beta1,beta2,gammax,sigma,logt){
             if (n<1){
               return(1)
             }else if (n<2){
               return(dnorm(logt[1],x%*%(beta1+beta2)+gammax+mu,sqrt(sigma)))
             }else{indx_list<-sapply(1:n,function(i){return(dnorm(logt[i],x[i,]%*%(beta1+beta2)+gammax+mu,sqrt(sigma)))})
             #indx<-prod(as.vector(indx_list))
             return(indx_list)}
           }
           
           f1_ind_cen<-function(n,mu,x,beta1,beta2,gammax,sigma,logt){
             if (n<1){
               return(1)
             }else if(n<2){
               return(1-pnorm(logt[1],x%*%(beta1+beta2)+gammax+mu,sqrt(sigma)))
             }else{    indx_list<-sapply(1:n,function(i){return(1-pnorm(logt[i],x[i,]%*%(beta1+beta2)+gammax+mu,sqrt(sigma)))})
             return(indx_list)}
             
           }
           
           f1_sum<-function(n,n_simu,mu,x_simu,beta1,beta2,gammax,sigma,log_t_sum){
             logt_mean<-apply(x_simu,2,mean)%*%(beta1+beta2)+gammax+mu
             logt_var<-sqrt(sigma/n)
             return(dnorm(log_t_sum,logt_mean,logt_var))
           } 
           
           
           t<-c()
           delta<-c()
           X<-c()
           tr<-c()
           
           for (i in 1:n_ipd){
             for (j in tr_list[[i]]){
               
               t<-c(t,t_list[[i]][[j]][[1]])
               delta<-c(delta,rep(1,n_list[[i]][[j]][[1]]))
               X<-rbind(X,X_list[[i]][[j]][[1]])
               tr<-c(tr,rep(j,n_list[[i]][[j]][[1]]))
               
               t<-c(t,t_list[[i]][[j]][[2]])
               delta<-c(delta,rep(1,n_list[[i]][[j]][[2]]))
               X<-rbind(X,X_list[[i]][[j]][[2]])
               tr<-c(tr,rep(j,n_list[[i]][[j]][[2]]))
               
             }
           }
           
           
           
           fre_res<-survreg(Surv(exp(t),delta)~X*tr,dist='lognormal')
           eta_fre<-fre_res$scale^2
           
           mu_mc<-mu_init
           sigma_mc<-sigma_init
           gamma_mc<-gamma_init
           
           beta1_mc<-beta1_init
           beta2_mc<-beta2_init
           
           
           mu_tab<-matrix(0,nrow=(burnin+nloop),ncol=ns)
           sigma_tab<-matrix(0,nrow=(burnin+nloop),ncol=ns)
           gamma_tab<-matrix(0,nrow=(burnin+nloop),ncol=nt)
           beta1_tab<-matrix(0,nrow=(burnin+nloop),ncol=n_cov)
           beta2_tab<-matrix(0,nrow=(burnin+nloop),ncol=n_cov)
           
           
           for (loop in 1:(burnin+nloop)){
             
             
             #mu:ipd
             
             mu_mc_tp<-sapply(1:n_ipd,function(i){
               mu_new<-rnorm(1,mu_mc[i],sd_pr)
               
               ind_ncen<-1
               if (tr_list[[i]][1]==1){
                 ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[1]][[1]],mu_new,X_list[[i]][[1]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[1]])/f1_ind_ncen(n_list[[i]][[1]][[1]],mu_mc[i],X_list[[i]][[1]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[1]]))
                 
                 for (j in tr_list[[i]][-1]){
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_new,X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                 }
               }else{
                 for (j in tr_list[[i]]){
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_new,X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                 }
               }
               
               ind_cen<-1
               if (tr_list[[i]][1]==1){
                 ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[1]][[2]],mu_new,X_list[[i]][[1]][[2]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[2]])/f1_ind_cen(n_list[[i]][[1]][[2]],mu_mc[i],X_list[[i]][[1]][[2]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[2]]))
                 
                 for (j in tr_list[[i]][-1]){
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_new,X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                 }
               }else{
                 for (j in tr_list[[i]]){
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_new,X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                 }
               }
               
               
               r<-ind_ncen*ind_cen
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,mu_new,mu_mc[i]))
             })
             
             mu_mc[1:n_ipd]<-mu_mc_tp
             
             #mu:agd
             
             mu_mc_tp<-sapply((n_ipd+1):ns,function(i){
               mu_new<-rnorm(1,mu_mc[i],sd_pr)
               
               r<-1
               
               if (tr_list[[i]][1]==1){
                 MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_new,X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]])
                 MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]])
                 r<-r*MH_new/MH_mc
                 
                 for (j in tr_list[[i]][-1]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   r<-r*MH_new/MH_mc
                 }
               }else{
                 for (j in tr_list[[i]]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   r<-r*MH_new/MH_mc
                 }
               }
               
               
               
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,mu_new,mu_mc[i]))
             })
             
             mu_mc[(n_ipd+1):ns]<-mu_mc_tp
             
             #sigma:ipd
             
             sigma_mc_tp<-sapply(1:n_ipd,function(i){
               sigma_new<-rnorm(1,eta_fre,0.1)
               
               ind_ncen<-1
               if (tr_list[[i]][1]==1){
                 ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[1]][[1]],mu_mc[i],X_list[[i]][[1]][[1]],beta1_mc,rep(0,n_cov),0,sigma_new,t_list[[i]][[1]][[1]])/f1_ind_ncen(n_list[[i]][[1]][[1]],mu_mc[i],X_list[[i]][[1]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[1]]))
                 
                 for (j in tr_list[[i]][-1]){
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                 }
               }else{
                 for (j in tr_list[[i]]){
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                 }
               }
               
               ind_cen<-1
               if (tr_list[[i]][1]==1){
                 ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[1]][[2]],mu_mc[i],X_list[[i]][[1]][[2]],beta1_mc,rep(0,n_cov),0,sigma_new,t_list[[i]][[1]][[2]])/f1_ind_cen(n_list[[i]][[1]][[2]],mu_mc[i],X_list[[i]][[1]][[2]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[2]]))
                 
                 for (j in  tr_list[[i]][-1]){
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                 }
               }else{
                 for (j in  tr_list[[i]]){
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                 }
               }
               
               
               r<-ind_ncen*ind_cen*dnorm(sigma_mc[i],eta_fre,0.1)/dnorm(sigma_new,eta_fre,0.1)
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,sigma_new,sigma_mc[i]))
             })
             
             sigma_mc[1:n_ipd]<-sigma_mc_tp
             
             
             #sigma:agd
             
             sigma_mc_tp<-sapply((n_ipd+1):ns,function(i){
               sigma_new<-rnorm(1,eta_fre,0.1)
               
               r<-1
               
               
               if (tr_list[[i]][1]==1){
                 MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,sigma_new,t_list[[i]][[1]])
                 MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]])
                 r<-r*MH_new/MH_mc
                 
                 for (j in  tr_list[[i]][-1]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   r<-r*MH_new/MH_mc
                 }
               }else{
                 for (j in  tr_list[[i]]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   r<-r*MH_new/MH_mc
                 }
               }
               
               
               
               
               r<-ifelse(is.na(r),0,min(r,1))*dnorm(sigma_mc[i],eta_fre,0.1)/dnorm(sigma_new,eta_fre,0.1)
               return(ifelse(runif(1)<r,sigma_new,sigma_mc[i]))
             })
             
             sigma_mc[(n_ipd+1):ns]<-sigma_mc_tp
             
             
             #gamma
             
             gamma_mc_tp<-sapply(2:nt,function(j){
               gamma_new<-rnorm(1,gamma_mc[j],0.01)
               
               ipd_s<-st_list[[j]][which(st_list[[j]]<=n_ipd)]
               agd_s<-st_list[[j]][which(st_list[[j]]>n_ipd)]
               
               ipd_r<-1
               if (length(ipd_s)>0){
                 ind_ncen<-1
                 ind_cen<-1
                 for (i in ipd_s){
                   jj<-j
                   
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[jj]][[1]],mu_mc[i],X_list[[i]][[jj]][[1]],beta1_mc,beta2_mc,gamma_new,sigma_mc[i],t_list[[i]][[jj]][[1]])/f1_ind_ncen(n_list[[i]][[jj]][[1]],mu_mc[i],X_list[[i]][[jj]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[jj]][[1]]))
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[jj]][[2]],mu_mc[i],X_list[[i]][[jj]][[2]],beta1_mc,beta2_mc,gamma_new,sigma_mc[i],t_list[[i]][[jj]][[2]])/f1_ind_cen(n_list[[i]][[jj]][[2]],mu_mc[1],X_list[[i]][[jj]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[jj]][[2]]))
                   ipd_r<-ipd_r*ind_ncen*ind_cen
                 }
               }
               
               agd_r<-1
               if (length(agd_s)>0){
                 
                 for (i in agd_s){
                   jj<-j
                   
                   MH_new<-f1_sum(n_list[[i]][[jj]],n_simu_list[[i]][[jj]],mu_mc[i],X_list[[i]][[jj]],beta1_mc,beta2_mc,gamma_new,sigma_mc[i],t_list[[i]][[jj]])
                   MH_mc<-f1_sum(n_list[[i]][[jj]],n_simu_list[[i]][[jj]],mu_mc[i],X_list[[i]][[jj]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[jj]])
                   agd_r<-agd_r*MH_new/MH_mc
                 }
               }
               
               r<-ipd_r*agd_r
               r<-ifelse(is.na(r),0,min(ipd_r*agd_r,1))
               return(ifelse(runif(1)<r,gamma_new,gamma_mc[j]))
             })
             
             gamma_mc[2:nt]<-gamma_mc_tp
             
             
             
             #beta1_mc
             
             
             for (l in 1:length(beta1_mc)){
               beta1_new<-beta1_mc
               beta1_new[l]<-rnorm(1,beta1_mc[l],0.1*sd_pr)
               
               beta1_ipd_tp<-sapply(1:n_ipd,function(i){
                 ind_ncen<-1
                 if (tr_list[[i]][1]==1){
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[1]][[1]],mu_mc[i],X_list[[i]][[1]][[1]],beta1_new,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[1]])/f1_ind_ncen(n_list[[i]][[1]][[1]],mu_mc[i],X_list[[i]][[1]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[1]]))
                   
                   for (j in tr_list[[i]][-1]){
                     ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                   }
                 }
                 
                 ind_cen<-1
                 if (tr_list[[i]][1]==1){
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[1]][[2]],mu_mc[i],X_list[[i]][[1]][[2]],beta1_new,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[2]])/f1_ind_cen(n_list[[i]][[1]][[2]],mu_mc[i],X_list[[i]][[1]][[2]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[2]]))
                   
                   for (j in tr_list[[i]][-1]){
                     ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                   }
                 }
                 return(ind_ncen*ind_cen)
               })
               
               beta1_agd_tp<-sapply((n_ipd+1):ns,function(i){
                 r<-1
                 
                 
                 if (tr_list[[i]][1]==1){
                   MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_new,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]])
                   MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]])
                   r<-r*MH_new/MH_mc
                   
                   for (j in tr_list[[i]][-1]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                     r<-r*MH_new/MH_mc
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                     r<-r*MH_new/MH_mc
                   }
                 }
                 
                 
                 
                 return(r)
                 
               })
               
               
               r<-prod(beta1_ipd_tp)*prod(beta1_agd_tp)
               r<-ifelse(is.na(r),0,min(r,1))
               beta1_mc[l]<-ifelse(runif(1)<r,beta1_new[l],beta1_mc[l])
               
               
             }
             
             #beta2_mc
             
             for (l in 1:length(beta2_mc)){
               beta2_new<-beta2_mc
               beta2_new[l]<-rnorm(1,beta2_mc[l],0.1*sd_pr)
               
               beta2_ipd_tp<-sapply(1:n_ipd,function(i){
                 ind_ncen<-1
                 if (tr_list[[i]][1]==1){
                   
                   
                   for (j in tr_list[[i]][-1]){
                     ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_new,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_new,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                   }
                 }
                 
                 ind_cen<-1
                 if (tr_list[[i]][1]==1){
                   
                   
                   for (j in tr_list[[i]][-1]){
                     ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_new,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[1],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_new,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[1],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                   }
                 }
                 return(ind_ncen*ind_cen)
               })
               
               beta2_agd_tp<-sapply((n_ipd+1):ns,function(i){
                 r<-1
                 for (j in 1:length(tr_list[[i]])){
                   
                   if (tr_list[[i]][1]==1){
                     
                     for (j in tr_list[[i]][-1]){
                       MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                       MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                       r<-r*MH_new/MH_mc
                     }
                   }else{
                     for (j in tr_list[[i]]){
                       MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                       MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                       r<-r*MH_new/MH_mc
                     }
                   }
                   
                 }
                 
                 return(r)
                 
               })
               
               r<-prod(beta2_ipd_tp)*prod(beta2_agd_tp)
               r<-ifelse(is.na(r),0,min(r,1))
               beta2_mc[l]<-ifelse(runif(1)<r,beta2_new[l],beta2_mc[l])
               
               
             }
             
             
             mu_tab[loop,]<-mu_mc
             sigma_tab[loop,]<-sigma_mc
             gamma_tab[loop,]<-gamma_mc
             beta1_tab[loop,]<-beta1_mc
             beta2_tab[loop,]<-beta2_mc
           }
           
           
           mu_est<-apply(mu_tab[(burnin+1):(burnin+nloop),],2,mean)
           sigma_est<-apply(sigma_tab[(burnin+1):(burnin+nloop),],2,mean)
           gamma_est<-apply(gamma_tab[(burnin+1):(burnin+nloop),],2,mean)
           beta1_est<-apply(beta1_tab[(burnin+1):(burnin+nloop),],2,mean)
           beta2_est<-apply(beta2_tab[(burnin+1):(burnin+nloop),],2,mean)
           
           return(list(mu_est,sigma_est,gamma_est,beta1_est,beta2_est,mu_tab[(burnin+1):(burnin+nloop),],sigma_tab[(burnin+1):(burnin+nloop),],gamma_tab[(burnin+1):(burnin+nloop),],beta1_tab[(burnin+1):(burnin+nloop),],beta2_tab[(burnin+1):(burnin+nloop),]))
           
         }, loglog={
           
           
           f1_ind_ncen<-function(n,mu,x,beta1,beta2,gammax,sigma,logt){
             
             if (n<1){
               return(1)
             }else if (n<2){
               return(dlogis(logt[1],x%*%(beta1+beta2)+gammax+mu,sigma))
             }else{
               indx_list<-sapply(1:n,function(i){return(dlogis(logt[i],x[i,]%*%(beta1+beta2)+gammax+mu,sigma))})
               #indx<-prod(as.vector(indx_list))
               return(indx_list)
             }
             
           }
           
           f1_ind_cen<-function(n,mu,x,beta1,beta2,gammax,sigma,logt){
             
             
             if (n<1){
               return(1)
             }else if (n<2){
               return(1-plogis(logt[1],x%*%(beta1+beta2)+gammax+mu,sigma))
             }else{
               indx_list<-sapply(1:n,function(i){return(1-plogis(logt[i],x[i,]%*%(beta1+beta2)+gammax+mu,sigma))})
               return(indx_list)
             }
             
           }
           
           f1_sum<-function(n,n_simu,mu,x_simu,beta1,beta2,gammax,sigma,log_t_sum){
             logt_mean<-apply(x_simu,2,mean)%*%(beta1+beta2)+gammax+mu
             logt_var<-sigma^2*pi/3/n
             return(dnorm(log_t_sum,logt_mean,logt_var))
           }
           
           
           t<-c()
           delta<-c()
           X<-c()
           tr<-c()
           
           for (i in 1:n_ipd){
             for (j in tr_list[[i]]){
               
               t<-c(t,t_list[[i]][[j]][[1]])
               delta<-c(delta,rep(1,n_list[[i]][[j]][[1]]))
               X<-rbind(X,X_list[[i]][[j]][[1]])
               tr<-c(tr,rep(j,n_list[[i]][[j]][[1]]))
               
               t<-c(t,t_list[[i]][[j]][[2]])
               delta<-c(delta,rep(1,n_list[[i]][[j]][[2]]))
               X<-rbind(X,X_list[[i]][[j]][[2]])
               tr<-c(tr,rep(j,n_list[[i]][[j]][[2]]))
               
             }
           }
           
           
           
           fre_res<-survreg(Surv(exp(t),delta)~X*tr,dist='loglogistic')
           eta_fre<-fre_res$scale
           
           
           mu_mc<-mu_init
           sigma_mc<-sigma_init
           gamma_mc<-gamma_init
           
           beta1_mc<-beta1_init
           beta2_mc<-beta2_init
           
           
           mu_tab<-matrix(0,nrow=(burnin+nloop),ncol=ns)
           sigma_tab<-matrix(0,nrow=(burnin+nloop),ncol=ns)
           gamma_tab<-matrix(0,nrow=(burnin+nloop),ncol=nt)
           beta1_tab<-matrix(0,nrow=(burnin+nloop),ncol=n_cov)
           beta2_tab<-matrix(0,nrow=(burnin+nloop),ncol=n_cov)
           
           
           for (loop in 1:(burnin+nloop)){
             
             
             #mu:ipd
             
             mu_mc_tp<-sapply(1:n_ipd,function(i){
               mu_new<-rnorm(1,mu_mc[i],sd_pr)
               
               ind_ncen<-1
               if (tr_list[[i]][1]==1){
                 ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[1]][[1]],mu_new,X_list[[i]][[1]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[1]])/f1_ind_ncen(n_list[[i]][[1]][[1]],mu_mc[i],X_list[[i]][[1]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[1]]))
                 
                 for (j in tr_list[[i]][-1]){
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_new,X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                 }
               }else{
                 for (j in tr_list[[i]]){
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_new,X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                 }
               }
               
               ind_cen<-1
               if (tr_list[[i]][1]==1){
                 ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[1]][[2]],mu_new,X_list[[i]][[1]][[2]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[2]])/f1_ind_cen(n_list[[i]][[1]][[2]],mu_mc[i],X_list[[i]][[1]][[2]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[2]]))
                 
                 for (j in tr_list[[i]][-1]){
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_new,X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                 }
               }else{
                 for (j in tr_list[[i]]){
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_new,X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                 }
               }
               
               
               r<-ind_ncen*ind_cen
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,mu_new,mu_mc[i]))
             })
             
             mu_mc[1:n_ipd]<-mu_mc_tp
             
             #mu:agd
             
             mu_mc_tp<-sapply((n_ipd+1):ns,function(i){
               mu_new<-rnorm(1,mu_mc[i],sd_pr)
               
               r<-1
               
               if (tr_list[[i]][1]==1){
                 MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_new,X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]])
                 MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]])
                 r<-r*MH_new/MH_mc
                 
                 for (j in tr_list[[i]][-1]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   r<-r*MH_new/MH_mc
                 }
               }else{
                 for (j in tr_list[[i]]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_new,X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   r<-r*MH_new/MH_mc
                 }
               }
               
               
               
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,mu_new,mu_mc[i]))
             })
             
             mu_mc[(n_ipd+1):ns]<-mu_mc_tp
             
             #sigma:ipd
             
             sigma_mc_tp<-sapply(1:n_ipd,function(i){
               sigma_new<-rnorm(1,eta_fre,0.1)
               
               ind_ncen<-1
               if (tr_list[[i]][1]==1){
                 ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[1]][[1]],mu_mc[i],X_list[[i]][[1]][[1]],beta1_mc,rep(0,n_cov),0,sigma_new,t_list[[i]][[1]][[1]])/f1_ind_ncen(n_list[[i]][[1]][[1]],mu_mc[i],X_list[[i]][[1]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[1]]))
                 
                 for (j in tr_list[[i]][-1]){
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                 }
               }else{
                 for (j in tr_list[[i]]){
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                 }
               }
               
               ind_cen<-1
               if (tr_list[[i]][1]==1){
                 ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[1]][[2]],mu_mc[i],X_list[[i]][[1]][[2]],beta1_mc,rep(0,n_cov),0,sigma_new,t_list[[i]][[1]][[2]])/f1_ind_cen(n_list[[i]][[1]][[2]],mu_mc[i],X_list[[i]][[1]][[2]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[2]]))
                 
                 for (j in  tr_list[[i]][-1]){
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                 }
               }else{
                 for (j in  tr_list[[i]]){
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                 }
               }
               
               
               r<-ind_ncen*ind_cen*dnorm(sigma_mc[i],eta_fre,0.1)/dnorm(sigma_new,eta_fre,0.1)
               r<-ifelse(is.na(r),0,min(r,1))
               return(ifelse(runif(1)<r,sigma_new,sigma_mc[i]))
             })
             
             sigma_mc[1:n_ipd]<-sigma_mc_tp
             
             
             #sigma:agd
             
             sigma_mc_tp<-sapply((n_ipd+1):ns,function(i){
               sigma_new<-rnorm(1,eta_fre,0.1)
               
               r<-1
               
               
               if (tr_list[[i]][1]==1){
                 MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,sigma_new,t_list[[i]][[1]])
                 MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]])
                 r<-r*MH_new/MH_mc
                 
                 for (j in  tr_list[[i]][-1]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   r<-r*MH_new/MH_mc
                 }
               }else{
                 for (j in  tr_list[[i]]){
                   MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_new,t_list[[i]][[j]])
                   MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                   r<-r*MH_new/MH_mc
                 }
               }
               
               
               
               
               r<-ifelse(is.na(r),0,min(r,1))*dnorm(sigma_mc[i],eta_fre,0.1)/dnorm(sigma_new,eta_fre,0.1)
               return(ifelse(runif(1)<r,sigma_new,sigma_mc[i]))
             })
             
             sigma_mc[(n_ipd+1):ns]<-sigma_mc_tp
             
             
             #gamma
             
             gamma_mc_tp<-sapply(2:nt,function(j){
               gamma_new<-rnorm(1,gamma_mc[j],0.01)
               
               ipd_s<-st_list[[j]][which(st_list[[j]]<=n_ipd)]
               agd_s<-st_list[[j]][which(st_list[[j]]>n_ipd)]
               
               ipd_r<-1
               if (length(ipd_s)>0){
                 ind_ncen<-1
                 ind_cen<-1
                 for (i in ipd_s){
                   jj<-j
                   
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[jj]][[1]],mu_mc[i],X_list[[i]][[jj]][[1]],beta1_mc,beta2_mc,gamma_new,sigma_mc[i],t_list[[i]][[jj]][[1]])/f1_ind_ncen(n_list[[i]][[jj]][[1]],mu_mc[i],X_list[[i]][[jj]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[jj]][[1]]))
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[jj]][[2]],mu_mc[i],X_list[[i]][[jj]][[2]],beta1_mc,beta2_mc,gamma_new,sigma_mc[i],t_list[[i]][[jj]][[2]])/f1_ind_cen(n_list[[i]][[jj]][[2]],mu_mc[1],X_list[[i]][[jj]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[jj]][[2]]))
                   ipd_r<-ipd_r*ind_ncen*ind_cen
                 }
               }
               
               agd_r<-1
               if (length(agd_s)>0){
                 
                 for (i in agd_s){
                   jj<-j
                   
                   MH_new<-f1_sum(n_list[[i]][[jj]],n_simu_list[[i]][[jj]],mu_mc[i],X_list[[i]][[jj]],beta1_mc,beta2_mc,gamma_new,sigma_mc[i],t_list[[i]][[jj]])
                   MH_mc<-f1_sum(n_list[[i]][[jj]],n_simu_list[[i]][[jj]],mu_mc[i],X_list[[i]][[jj]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[jj]])
                   agd_r<-agd_r*MH_new/MH_mc
                 }
               }
               
               r<-ipd_r*agd_r
               r<-ifelse(is.na(r),0,min(ipd_r*agd_r,1))
               return(ifelse(runif(1)<r,gamma_new,gamma_mc[j]))
             })
             
             gamma_mc[2:nt]<-gamma_mc_tp
             
             
             
             #beta1_mc
             
             
             for (l in 1:length(beta1_mc)){
               beta1_new<-beta1_mc
               beta1_new[l]<-rnorm(1,beta1_mc[l],0.1*sd_pr)
               
               beta1_ipd_tp<-sapply(1:n_ipd,function(i){
                 ind_ncen<-1
                 if (tr_list[[i]][1]==1){
                   ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[1]][[1]],mu_mc[i],X_list[[i]][[1]][[1]],beta1_new,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[1]])/f1_ind_ncen(n_list[[i]][[1]][[1]],mu_mc[i],X_list[[i]][[1]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[1]]))
                   
                   for (j in tr_list[[i]][-1]){
                     ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                   }
                 }
                 
                 ind_cen<-1
                 if (tr_list[[i]][1]==1){
                   ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[1]][[2]],mu_mc[i],X_list[[i]][[1]][[2]],beta1_new,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[2]])/f1_ind_cen(n_list[[i]][[1]][[2]],mu_mc[i],X_list[[i]][[1]][[2]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]][[2]]))
                   
                   for (j in tr_list[[i]][-1]){
                     ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                   }
                 }
                 return(ind_ncen*ind_cen)
               })
               
               beta1_agd_tp<-sapply((n_ipd+1):ns,function(i){
                 r<-1
                 
                 
                 if (tr_list[[i]][1]==1){
                   MH_new<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_new,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]])
                   MH_mc<-f1_sum(n_list[[i]][[1]],n_simu_list[[i]][[1]],mu_mc[i],X_list[[i]][[1]],beta1_mc,rep(0,n_cov),0,sigma_mc[i],t_list[[i]][[1]])
                   r<-r*MH_new/MH_mc
                   
                   for (j in tr_list[[i]][-1]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                     r<-r*MH_new/MH_mc
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                     MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                     r<-r*MH_new/MH_mc
                   }
                 }
                 
                 
                 
                 return(r)
                 
               })
               
               
               r<-prod(beta1_ipd_tp)*prod(beta1_agd_tp)
               r<-ifelse(is.na(r),0,min(r,1))
               beta1_mc[l]<-ifelse(runif(1)<r,beta1_new[l],beta1_mc[l])
               
               
             }
             
             #beta2_mc
             
             for (l in 1:length(beta2_mc)){
               beta2_new<-beta2_mc
               beta2_new[l]<-rnorm(1,beta2_mc[l],0.1*sd_pr)
               
               beta2_ipd_tp<-sapply(1:n_ipd,function(i){
                 ind_ncen<-1
                 if (tr_list[[i]][1]==1){
                   
                   
                   for (j in tr_list[[i]][-1]){
                     ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_new,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ind_ncen<-ind_ncen*prod(f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_new,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]])/f1_ind_ncen(n_list[[i]][[j]][[1]],mu_mc[i],X_list[[i]][[j]][[1]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[1]]))
                   }
                 }
                 
                 ind_cen<-1
                 if (tr_list[[i]][1]==1){
                   
                   
                   for (j in tr_list[[i]][-1]){
                     ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_new,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[1],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                   }
                 }else{
                   for (j in tr_list[[i]]){
                     ind_cen<-ind_cen*prod(f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[i],X_list[[i]][[j]][[2]],beta1_mc,beta2_new,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]])/f1_ind_cen(n_list[[i]][[j]][[2]],mu_mc[1],X_list[[i]][[j]][[2]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]][[2]]))
                   }
                 }
                 return(ind_ncen*ind_cen)
               })
               
               beta2_agd_tp<-sapply((n_ipd+1):ns,function(i){
                 r<-1
                 for (j in 1:length(tr_list[[i]])){
                   
                   if (tr_list[[i]][1]==1){
                     
                     for (j in tr_list[[i]][-1]){
                       MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                       MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                       r<-r*MH_new/MH_mc
                     }
                   }else{
                     for (j in tr_list[[i]]){
                       MH_new<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_new,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                       MH_mc<-f1_sum(n_list[[i]][[j]],n_simu_list[[i]][[j]],mu_mc[i],X_list[[i]][[j]],beta1_mc,beta2_mc,gamma_mc[j],sigma_mc[i],t_list[[i]][[j]])
                       r<-r*MH_new/MH_mc
                     }
                   }
                   
                 }
                 
                 return(r)
                 
               })
               
               r<-prod(beta2_ipd_tp)*prod(beta2_agd_tp)
               r<-ifelse(is.na(r),0,min(r,1))
               beta2_mc[l]<-ifelse(runif(1)<r,beta2_new[l],beta2_mc[l])
               
               
             }
             
             
             mu_tab[loop,]<-mu_mc
             sigma_tab[loop,]<-sigma_mc
             gamma_tab[loop,]<-gamma_mc
             beta1_tab[loop,]<-beta1_mc
             beta2_tab[loop,]<-beta2_mc
           }
             

           
           mu_est<-apply(mu_tab[(burnin+1):(burnin+nloop),],2,mean)
           sigma_est<-apply(sigma_tab[(burnin+1):(burnin+nloop),],2,mean)
           gamma_est<-apply(gamma_tab[(burnin+1):(burnin+nloop),],2,mean)
           beta1_est<-apply(beta1_tab[(burnin+1):(burnin+nloop),],2,mean)
           beta2_est<-apply(beta2_tab[(burnin+1):(burnin+nloop),],2,mean)
           
           return(list(mu_est,sigma_est,gamma_est,beta1_est,beta2_est,mu_tab[(burnin+1):(burnin+nloop),],sigma_tab[(burnin+1):(burnin+nloop),],gamma_tab[(burnin+1):(burnin+nloop),],beta1_tab[(burnin+1):(burnin+nloop),],beta2_tab[(burnin+1):(burnin+nloop),]))
         })
  
}
