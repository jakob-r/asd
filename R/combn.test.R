combn.test <-
function(stage1,stage2,weight=0.5,method="invnorm"){

 # valid test type
 meths <- c("invnorm", "fisher")
 imeths <- as.integer(match(method, meths, -1))
 if (imeths < 1) {
  stop("Invalid combination test type")
 } # end if

 # extract pvalues from dunnett test object
 p1 <- stage1$pvalues
 p2 <- stage2$pvalues

 # intersection hypotheses indicator
 hyp.comb <- stage1$hyp.comb

 # number of treatments is length of Z vector
 treats <- length(p1)

 # weights
 if (method=="invnorm"){ 
  w1 <- sqrt(weight)
  w2 <- sqrt(1-w1*w1)
 } else {
  w1 <- 1; w2 <- 1
 } # end if

 # structure zscores same as p1
 zscores <- p1
 
 # identify elements of zscores
 for (i in 1:treats){
  for (j in 1:length(p1[[i]])){
  
   # replace with missing value indicator
   zscores[[i]][j] <- NA
   
   # inverse normal combination
   if (method=="invnorm"){
    zscores[[i]][j] <- w1*qnorm(1-p1[[i]][j])+w2*qnorm(1-p2[[i]][j])
   } # end if

   # fisher combination
   if (method=="fisher"){
    zscores[[i]][j] <- qnorm(pchisq(-2*log(p1[[i]][j]*p2[[i]][j]),4))
   } # end if

  } # end j
 } # end i

 # results of combination test
 list(method=method,zscores=zscores,hyp.comb=hyp.comb,weights=c(w1,w2))

} # end of function combn.test

