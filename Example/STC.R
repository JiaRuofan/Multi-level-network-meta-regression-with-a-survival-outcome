AB_IPD <- cbind(x_AB,t1_censored,delta1,c(rep(0,n_AB/2),rep(1,n_AB/2)))
colnames(AB_IPD)=c('X1','X2','X3','X4','X5',t,delta,trt)

AC_AGD <- apply(x_AC,2,mean)
colnames(AC_AGD)=c('mean_X1','mean_X2','mean_X3','mean_X4','mean_X5')

# fit regression of outcome on the baseline characteristics and treatment
# effect modifiers are centered at the mean AC values

outcome.fit <- coxph(Surv(t, delta)~X1+X2+X3+X4+X5+trt*I(X1-AC_AGD$mean_X1)+trt*I(X2-AC_AGD$mean_X2)+trt*I(X3-AC_AGD$mean_X3)+trt*I(X3-AC_AGD$mean_X3)+trt*I(X3-AC_AGD$mean_X3),
                     data=AB_IPD)
