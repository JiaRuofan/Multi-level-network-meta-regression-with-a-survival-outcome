AB_IPD <- read.csv("AB_IPD.csv") # load AB patient-level data
AC_AGD <- read.csv("AC_AGD.csv") # load AC aggregate-level data

# fit regression of outcome on the baseline characteristics and treatment
# effect modifiers are centered at the mean AC values

outcome.fit <- coxph(Surv(t, delta)~X1+X2+X3+X4+X5+trt*I(X1-AC_AGD$mean_X1)+trt*I(X2-AC_AGD$mean_X2)+trt*I(X3-AC_AGD$mean_X3)+trt*I(X3-AC_AGD$mean_X3)+trt*I(X3-AC_AGD$mean_X3),
                     data=AB_IPD)
