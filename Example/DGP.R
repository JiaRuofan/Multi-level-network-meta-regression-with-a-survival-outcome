
ns<-2 #Number of studies
nt<-3 #Number of treatments
n_cov<-5 #Number of covariates

n_AB<-200 #sample size
n_AC<-200

beta1<-c(0.1,0.1,-0.1,0.05,-0.25)
#beta2<-c(-0.02,-0.02,-0.1)
beta2<-c(-0.2,0.2,-0.2,-0.2,-0.1)

gammaB<--1.2
gammaC<--0.5

eta1<-0.8
eta2<-1.2

mu1<-6.2
mu2<-5.8

shape1<-eta1
shape2<-eta2

n1_A<-n_AB/2;n1_B<-n_AB/2;n2_A<-n_AC/2;n2_C<-n_AC/2

sigma=diag(1,3)
sigma[1,2]<-sigma[2,1]<-sigma[3,2]<-sigma[2,3]<-0.5
sigma[1,3]<-sigma[3,1]<-0.1
x_AB<-cbind(mvrnorm(n_AB,mu=rep(0,3),sigma),rgamma(n_AB,4,2),rbinom(n_AB,1,0.2))
sigma=diag(1,3)
sigma[1,2]<-sigma[2,1]<-sigma[3,2]<-sigma[2,3]<-0.5
sigma[1,3]<-sigma[3,1]<-0.1
x_AC<-cbind(mvrnorm(n_AC,mu=rep(0.5,3),sigma),rgamma(n_AC,6,2),rbinom(n_AC,1,0.7))
   
  
  
scale1<-1/c(mu1^(1/eta1)*exp((x_AB[1:(n_AB/2),]%*%(beta1))/eta1),mu1^(1/eta1)*exp((x_AB[(n_AB/2+1):(n_AB),]%*%(beta1+beta2)+gammaB)/eta1))
scale2<-1/c(mu2^(1/eta2)*exp((x_AC[1:(n_AC/2),]%*%(beta1))/eta2),mu2^(1/eta2)*exp((x_AC[(n_AC/2+1):(n_AC),]%*%(beta1+beta2)+gammaC)/eta2))
    
t1<-c()
  for (i in 1:n_AB){
    
    t1[i]<-rweibull(1,shape=shape1,scale=scale1[i])
}
  
  
t2<-c()
for (i in 1:n_AC){
    
    t2[i]<-rweibull(1,shape=shape2,scale=scale2[i])
}
  
t2_A_sum<-mean(t2[1:n2_A])
t2_C_sum<-mean(t2[(n2_A+1):(n2_A+n2_C)])
  
    
    
t1_censored<-c()
delta1<-c()
  
for (i in 1:n_AB){
    if (t1[i]>0.6){
      t1_censored[i]<-0.6
      delta1[i]<-0
    }else{
      t1_censored[i]<-t1[i]
      delta1[i]<-1
    }
}
  
#covariates  
x1_A<-x_AB[1:(n_AB/2),]
x1_B<-x_AB[(n_AB/2+1):n_AB,]
  
x2_A<-x_AC[1:(n_AC/2),]
x2_C<-x_AC[(n_AC/2+1):n_AC,]
    
  
#t 
t1_A<-t1_censored[1:n1_A]
t1_B<-t1_censored[(n1_A+1):(n1_A+n1_B)]

#censoring
delta1_A<-delta1[1:n1_A]
delta1_B<-delta1[(n1_A+1):(n1_A+n1_B)]
