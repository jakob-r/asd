gasd.sim <-
function(z1=c(0,0,0),z2=c(0,0,0),zearly=c(0,0,0),
                  v1=c(1,1,1,1),v2=c(1,1,1,1),vearly=c(1,1,1,1),
                  corr=c(0,0,0,0),weight=0.5,nsim=1000,
                  seed=12345678,select=0,epsilon=1,thresh=1,
                  level=0.025,ptest=seq(1:length(z1)),
                  fu=FALSE,method="invnorm"){

 # validate inputs

 # z1, z2 and early
 ntreat <- length(z1)
 if(ntreat>8){
  stop("Maximum of 8 test treatments allowed")
 } # end if
 if(length(z2)!=ntreat){
  stop("z1 and z2 vectors must be same length")
 } # end if
 if(length(zearly)!=ntreat){
  stop("z1 and zearly vectors must be same length")
 } # end if

 # v1, v2 and vearly
 vntreat <- ntreat+1
 if(length(v1)!=vntreat){
  stop("v1 must be vector length(z1)+1")
 } # end if
 if(min(v1)<=0){
   stop("All elements of v1 must be >= 0")
 } # end if
 if(length(v2)!=vntreat){
  stop("v2 must be vector length(z1)+1")
 } # end if
 if(min(v2)<=0){
   stop("All elements of v2 must be >= 0")
 } # end if
 if(length(vearly)!=vntreat){
  stop("vearly must be vector length(z1)+1")
 } # end if
 if(min(vearly)<=0){
   stop("All elements of vearly must be >= 0")
 } # end if

 # correlation
 if (length(corr)==vntreat){
  corext <- c(abs(max(corr)),abs(min(corr)))
  if (max(corext)>1){
   stop("Correlation must be between -1 and 1") 
  } # end if
 } # end if
 if(length(corr)!=vntreat){
   stop("Correlation must be vector length(z1)+1") 
 } # end if

 # seed
 set.seed(abs(round(seed,0)))

 # number of simulations
 nsim <- abs(round(nsim,0))
 if (nsim>10000000){
  stop("Maximum of 10,000,000 simulations allowed")
 } # end if

 # selection rule
 select <- abs(as.integer(select))
 if (select>6){
  stop("Selection rule options: 0,1,2,3,4,5,6")
 } # end if

 # level
 if (level>=1 | level<=0){
   stop("Level must be between 0 and 1") 
 } # end if
 
 # power tests
 if (length(ptest)>ntreat){
  stop("Invalid ptest vector: must be <= ntreat")
 } # end if
 ptest.check <- sum(is.element(ptest,seq(1:length(z1))))
 if(ptest.check<length(ptest)){
  stop("Invalid ptest vector: must be element of seq(1:length(z1))")
 } # end if

 # followup
 if(is.logical(fu)){
   follow <- fu
  } else {
   follow <- FALSE
 } # end if

 # method
 meth.options <- c("invnorm", "fisher")
 imeth <- as.integer(match(method,meth.options,-1))
 if (imeth<1) {
  stop("Unknown method: current options invnorm or fisher")
 } # end if


 # set-up simulations

 # structures to store results
 sim.reject <- 0
 select_total <- NULL
 reject_total <- rep(0,times=ntreat)
 count_total <- rep(0,times=ntreat)

 # construction
 lambda.1 <- sqrt(v1[1])/sqrt(v1[1]+v1[2:(ntreat+1)])
 lambda.2 <- sqrt(v2[1])/sqrt(v2[1]+v2[2:(ntreat+1)])
 lambda.early <- sqrt(vearly[1])/sqrt(vearly[1]+vearly[2:(ntreat+1)])
 lambda.1 <- matrix(lambda.1,ncol=1,nrow=ntreat)
 lambda.2 <- matrix(lambda.2,ncol=1,nrow=ntreat)
 lambda.early <- matrix(lambda.early,ncol=1,nrow=ntreat)

 # var(Z1,Z1)
 varm.1 <- lambda.1%*%t(lambda.1)+diag(c(1-(lambda.1)^2))

 # var(Z2,Z2)
 varm.2 <- lambda.2%*%t(lambda.2)+diag(c(1-(lambda.2)^2))

 # var(Z1*,Z1*)
 varm.early <- lambda.early%*%t(lambda.early)+diag(c(1-(lambda.early)^2))

 # cov(Z1,Z1*)
 mat.1 <- corr[1]*(lambda.early%*%t(lambda.1))
 fact.1 <- sqrt(v1[2:(ntreat+1)])/sqrt(v1[1])
 fact.early <- sqrt(vearly[2:(ntreat+1)])/sqrt(vearly[1])
 lambda.product <- lambda.1*lambda.early
 mat.2 <- diag(c(corr[2:(ntreat+1)]*fact.1*fact.early*lambda.product))
 cov.mat <- mat.1+mat.2 

 # full covariance matrix
 z.all <- c(z1,zearly,z2)
 munit <- length(z1)
 z.covm <- matrix(0,nrow=length(z.all),ncol=length(z.all))
 z.covm[1:munit,1:munit] <- varm.1
 z.covm[(munit+1):(2*munit),(munit+1):(2*munit)] <- varm.early
 z.covm[(2*munit+1):(3*munit),(2*munit+1):(3*munit)] <- varm.2
 z.covm[(munit+1):(2*munit),1:munit] <- cov.mat
 z.covm[1:munit,(munit+1):(2*munit)] <- t(cov.mat)


 # start simulations
 for (i in 1:nsim){

  # generate single randomization
  randi <- rmvnorm(n=1,mean=z.all,sigma=z.covm)

  # test statistics
  # early outcome
  Etest.stat <- c(as.numeric(randi[(munit+1):(2*munit)]))
  # stage 1
  Itest.stat <- c(as.numeric(randi[1:munit]))
  # stage 2
  Ftest.stat <- c(as.numeric(randi[(2*munit+1):length(randi)]))

  # treatment selection
  interim.select <- select.rule(x=Etest.stat,type=select,epsilon=epsilon,thresh=thresh)
  if (i==1){
   select_total <- interim.select$select} else {
   select_total <- select_total+interim.select$select
  } # end if

  # first stage p-values
  if (follow==TRUE){
   dunnett1 <- dunnett.test(Itest.stat)
  } else { 
   Itest.stat <- Itest.stat+interim.select$z
   dunnett1 <- dunnett.test(Itest.stat)
  } # end if

  # if treatments selected
  if (sum(interim.select$select)>0){
  
    # count selections
    count_i <- rep(0,times=ntreat)
    count_i[sum(interim.select$select)] <- 1
    count_total <- count_total+count_i

    # second stage p-values
    dunnett2 <- dunnett.test(Ftest.stat,select=interim.select$select)

    # combination test
    comb.test <- combn.test(dunnett1,dunnett2,weight=weight,method=method)

    # hypothesis testing
    reject_test <- hyp.test(comb.test,level=level,full.hyp=FALSE)
    if (i==1){
     reject_total <- reject_test$reject} else {
     reject_total <- reject_total+reject_test$reject
    } # end if

    # count if one or more of select treatments are rejected
    test.stat <- sum(reject_test$reject[ptest])
    if (test.stat>=1){
     sim.reject <- sim.reject+1
    } # end test of rejection

   } else {

    if (follow==TRUE){
     
     # follow-up if none selected

     # hypothesis testing
     reject_test <- hyp.test(dunnett1,level=level,full.hyp=FALSE)
     if (i==1){
      reject_total <- reject_test$reject} else {
      reject_total <- reject_total+reject_test$reject
     } # end if

     # count if one or more of select treatments are rejected
     test.stat <- sum(reject_test$reject[ptest])
     if (test.stat>=1){
      sim.reject <- sim.reject+1
     } # end test of rejection

    } else {

     # do not follow-up if none selected

     reject_total <- reject_total+rep(0,times=ntreat)

    } # end if

   } # end treatments selected


 } # end i simulations


 # results
 reject_total <- matrix(reject_total,ncol=ntreat,nrow=1)
 colnames(reject_total) <- colnames(reject_test$reject)
 rownames(reject_total) <- as.character("n")
 select_total <- matrix(select_total,ncol=ntreat,nrow=1)
 colnames(select_total) <- as.character(1:ntreat)
 rownames(select_total) <- as.character("n")
 count_total <- matrix(count_total,ncol=ntreat,nrow=1)
 colnames(count_total) <- as.character(1:ntreat)
 rownames(count_total) <- as.character("n")
 sim.reject <- matrix(sim.reject,ncol=1,nrow=1)
 colnames(sim.reject) <- "Total"
 rownames(sim.reject) <- as.character("n")
 list(count.total=count_total,select.total=select_total,
     reject.total=reject_total,sim.reject=sim.reject)


} # end of function gasd.sim

