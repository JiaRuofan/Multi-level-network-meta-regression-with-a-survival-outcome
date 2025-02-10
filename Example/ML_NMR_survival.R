################Function for simulating psedo-covariates
 
 cov_simu5<-function(x_dat,cor_refer,n_simu){
 n_cov<-dim(x_dat)[2]
 x1_mean<-mean(x_dat[,1])
 x1_var<-var(x_dat[,1])
 x2_mean<-mean(x_dat[,2])
 x2_var<-var(x_dat[,2])
 x3_mean<-mean(x_dat[,3])
 x3_var<-var(x_dat[,3])
 x4_mean<-mean(x_dat[,4])
 x4_var<-var(x_dat[,4])
 x5_p<-x5_mean<-mean(x_dat[,5])
 x4_scale<-x4_var/x4_mean
 x4_shape<-x4_mean/x4_scale
 ## Initialization and parameters
 P <- cor_refer
 d <- nrow(P)
 n_simu <- n_simu
 # Dimension
 # Number of samples
 ## Simulation (non-vectorized version)
 A <- t(chol(P))
 U <- matrix(nrow = n_simu, ncol = d)
 for (i in 1:n_simu){
 Z
 <- rnorm(d)
 X
 <- A%*%Z
 U[i, ] <- pnorm(X)
 }
 ## Simulation (compact vectorized version)
 U <- pnorm(matrix(rnorm(n_simu*d), ncol = d) %*% chol(P))
 X_simu<-matrix(0,nrow=n_simu,ncol=n_cov)
 for (i in 1:n_simu){
 X_simu[i,1]<-qnorm(U[i,1],x1_mean,sqrt(x1_var))
 X_simu[i,2]<-qnorm(U[i,2],x2_mean,sqrt(x2_var))
 X_simu[i,3]<-qnorm(U[i,3],x3_mean,sqrt(x3_var))
 X_simu[i,4]<-qgamma(U[i,4],x4_shape,scale=x4_scale)
 X_simu[i,5]<-ifelse(U[i,5]<x5_p,1,0)
 }
 return(X_simu)
 }


 #Construct parameters for mlnmr_survival

  ns=2
 n_cov=5
 nt=3
 n_ipd=1
 n_agd=1
 type=’weibull_median’
 n_simu=10000 #Sample size of simulated X
 burnin<-20000
 nloop<-10000

  #Simulate X for AgD
 n2_simu_A<-n_simu;n2_simu_C<-n_simu;
 X2_simu_A<-cov_simu5(x2_A,cor(x_AB),n2_simu_A);X2_simu_C<-cov_simu5(x2_C,cor(x_AB),n2_simu_C)
 n_list<-list()
 n_list[[1]]<-list()
 n_list[[1]][[1]]<-n1_A

 n_list[[1]][[2]]<-n1_B
 n_list[[2]]<-list()
 n_list[[2]][[1]]<-n2_A
 n_list[[2]][[3]]<-n2_C
 X_list<-list()
 X_list[[1]]<-list()
 X_list[[1]][[1]]<-x1_A
 X_list[[1]][[2]]<-x1_B
 X_list[[2]]<-list()
 X_list[[2]][[1]]<-X2_simu_A
 X_list[[2]][[3]]<-X2_simu_C
 n_simu_list<-list()
 n_simu_list[[2]]<-list()
 n_simu_list[[2]][[1]]<-n2_simu_A
 n_simu_list[[2]][[3]]<-n2_simu_C
 t_list<-list()
 t_list[[1]]<-list()
 t_list[[2]]<-list()
 t_list[[1]][[1]]<-t1_A
 t_list[[1]][[2]]<-t1_B
 t_list[[2]][[1]]<-mean(t2[1:n2_A])
 t_list[[2]][[3]]<-mean(t2[(n2_A+1):(n2_A+n2_C)])
 delta_list<-list()
 delta_list[[1]]<-list()
 delta_list[[1]][[1]]<-delta1_A
 delta_list[[1]][[2]]<-delta1_B
 tr_list<-list()
 tr_list[[1]]<-c(1,2)
 tr_list[[2]]<-c(1,3)
 st_list<-list()
 st_list[[1]]<-c(1,2,3,4)
 st_list[[2]]<-c(1)
 st_list[[3]]<-c(2)

mlnmr_survival(ns,n_ipd,n_agd,n_cov,nt,n_list,t_list,X_list,n_simu_list,delta_list,var_list,tr_list,st_list,burnin=20000,nloop=10000,sd_pr=0.1,type)

#Construct parameters for mlnmr_survival_trans


 ns=2
 n_cov=5
 nt=3
 n_ipd=1
 n_agd=1

 type=’logn’
 n_simu=10000 #sample size of simulated X
 burnin<-20000
 nloop<-10000
 #take log
 t1_censored<-log(t1_censored)
 #separate censored samples and non-censored samples
 t1_A_ncensored<-t1_censored[1:(n_AB/2)][which(delta1[1:(n_AB/2)]==1)]
 t1_A_censored<-t1_censored[1:(n_AB/2)][which(delta1[1:(n_AB/2)]==0)]
 t1_B_ncensored<-t1_censored[(n_AB/2+1):n_AB][which(delta1[(n_AB/2+1):n_AB]==1)]
 t1_B_censored<-t1_censored[(n_AB/2+1):n_AB][which(delta1[(n_AB/2+1):n_AB]==0)]
 x1_A_ncen<-x_AB[which(delta1[1:(n_AB/2)]==1),]
 x1_A_cen<-x_AB[which(delta1[1:(n_AB/2)]==0),]
 x1_B_ncen<-x_AB[(n_AB/2+1):n_AB,][which(delta1[(n_AB/2+1):n_AB]==1),]
 x1_B_cen<-x_AB[(n_AB/2+1):n_AB,][which(delta1[(n_AB/2+1):n_AB]==0),]
 x2_A<-x_AC[1:(n_AC/2),]
 x2_C<-x_AC[(n_AC/2+1):n_AC,]
 n1_A_ncen<-length(t1_A_ncensored);n1_A_cen<-length(t1_A_censored)
 n1_B_ncen<-length(t1_B_ncensored);n1_B_cen<-length(t1_B_censored)
 #Simulate X for AgD
 n2_simu_A<-n_simu;n2_simu_C<-n_simu;
 X2_simu_A<-cov_simu5(x2_A,cor(x_AB),n2_simu_A);X2_simu_C<-cov_simu5(x2_C,cor(x_AB),n
 n_list<-list()
 n_list[[1]]<-list()
 n_list[[1]][[1]]<-list()
 n_list[[1]][[2]]<-list()
 n_list[[1]][[1]][[1]]<-n1_A_ncen
 n_list[[1]][[1]][[2]]<-n1_A_cen
 n_list[[1]][[2]][[1]]<-n1_B_ncen
 n_list[[1]][[2]][[2]]<-n1_B_cen
 n_list[[2]]<-list()
 n_list[[2]][[1]]<-n2_A
 n_list[[2]][[3]]<-n2_C
 X_list<-list()
 X_list[[1]]<-list()
 X_list[[1]][[1]]<-list()
 X_list[[1]][[2]]<-list()
 X_list[[1]][[1]][[1]]<-x1_A_ncen
 X_list[[1]][[1]][[2]]<-x1_A_cen
 X_list[[1]][[2]][[1]]<-x1_B_ncen
 X_list[[1]][[2]][[2]]<-x1_B_cen
 X_list[[2]]<-list()
 X_list[[2]][[1]]<-X2_simu_A

 X_list[[2]][[3]]<-X2_simu_C
 n_simu_list<-list()
 n_simu_list[[2]]<-list()
 n_simu_list[[2]][[1]]<-n2_simu_A
 n_simu_list[[2]][[3]]<-n2_simu_C
 t_list<-list()
 t_list[[1]]<-list()
 t_list[[2]]<-list()
 t_list[[1]][[1]]<-list()
 t_list[[1]][[2]]<-list()
 t_list[[1]][[1]][[1]]<-t1_A_ncensored
 t_list[[1]][[1]][[2]]<-t1_A_censored
 t_list[[1]][[2]][[1]]<-t1_B_ncensored
 t_list[[1]][[2]][[2]]<-t1_B_censored
 t_list[[2]][[1]]<-mean(log(t2[1:n2_A]))
 t_list[[2]][[3]]<-mean(log(t2[(n2_A+1):(n2_A+n2_C)]))
 tr_list<-list()
 tr_list[[1]]<-c(1,2)
 tr_list[[2]]<-c(1,3)
 st_list<-list()
 st_list[[1]]<-c(1,2,3,4)
 st_list[[2]]<-c(1)
 st_list[[3]]<-c(2)

 mlnmr_survival_trans(ns,n_ipd,n_agd,n_cov,nt,n_list,t_list,X_list,n_simu_list,tr_list,st_list,burnin=20000,nloop=10000,sd_pr=0.1,type)
