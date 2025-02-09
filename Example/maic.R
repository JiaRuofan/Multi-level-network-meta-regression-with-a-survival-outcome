 #install.packages(’maic’)
 library(maic)

 data_AB<-cbind(rep(0,n1_A),rep(0,n1_B),x_AB)
 colnames(data_AB)<-c(’tr’,colnames(x_AB))

 target <- c("x1_mean" = mean(x_AC[,1]),
 "x2_mean" = mean(x_AC[,2]),
 "x3_mean" = mean(x_AC[,3]),
 "x4_mean" = mean(x_AC[,4]),
 "x5_mean" = mean(x_AC[,5]))

 dict <- data.frame(
 "match.id" =
 c("x1", "x2",
 "x3", "x4", "x5"),"target.variable" =
 c("x1_mean", "x2_mean","x3_mean","x4_mean","x5_mean"),
 "index.variable" =
 c("x1", "x2",
 "x3",
 "x4",
 "x5"),
 "match.type" =
 c("mean", "mean", "mean", "mean", "mean"),
 stringsAsFactors = FALSE)

 ipmat <- createMAICInput(
 index = x_AB,
 target = target,
 dictionary = dict,
 matching.variables =
 c("x1", "x2",
 "x3", "x4", "x5"))

 wts <- maicWeight(ipmat) #weights of samples of AB in AC
