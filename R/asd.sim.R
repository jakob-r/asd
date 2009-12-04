asd.sim <-
function(nsamp=c(32,32),early=c(0,0,0),final=c(0,0,0),nsim=1000,corr=0,
                     seed=12345678,select=0,epsilon=1,thresh=1,
                     level=0.025,ptest=seq(1:length(early)),reall=FALSE,
                     fu=FALSE,estim=FALSE,method="invnorm"){

 # validate inputs

 # sample size
 if(length(nsamp)!=2){
  stop("Sample size for stages 1 and 2 required: input vector length=2")
 } # end if
 if(length(nsamp)==2){
  nsamp <- abs(round(as.numeric(nsamp),0))
 } # end if

 # correlation
 if (length(corr)==1){
  if (abs(corr)>1){
   stop("Correlation must be between -1 and 1") 
  } # end if
 } # end if
 if(length(corr)>1){
   stop("Correlation must be vector length=1") 
 } # end if

 # seed
 set.seed(abs(round(seed,0)))

 # early and final
 ntreat <- length(early)
 if(ntreat>8){
  stop("Maximum of 8 test treatments allowed")
 } # end if
 if(length(final)!=ntreat){
  stop("Early and final vectors must be same length")
 } # end if

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
 ptest.check <- sum(is.element(ptest,seq(1:length(early))))
 if(ptest.check<length(ptest)){
  stop("Invalid ptest vector: must be element of seq(1:length(early))")
 } # end if

 # estimation
 if(is.logical(estim)){
   estim <- estim
  } else {
   estim <- FALSE
 } # end if

 # reallocation and followup
 if(is.logical(reall) && (select==4 | select==6)){
   real <- reall
  } else {
   real <- FALSE
 } # end if 
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

 # sample sizes
 if(real==TRUE){
  nsamp[1] <- floor(nsamp[1]/(ntreat+1))
  nsamp_stage2 <- nsamp[2]
  # to ensure weights are correct for comb.test
  nsamp[2] <- floor(nsamp[2]/(ntreat+1))
 } # end if


 # structures to store results
 sim.reject <- 0
 select_total <- NULL
 reject_total <- rep(0,times=ntreat)
 count_total <- rep(0,times=ntreat)
 if (estim==TRUE){
  sum1treat_1 <- rep(0,times=ntreat); sum1treat_2 <- rep(0,times=ntreat)
  sum2treat_1 <- rep(0,times=ntreat); sum2treat_2 <- rep(0,times=ntreat)
  sumN_1 <- rep(0,times=ntreat); sumN_2 <- rep(0,times=ntreat)
  est_mean <- matrix(NA,ncol=ntreat,nrow=2)
  est_N <- matrix(NA,ncol=ntreat,nrow=2)
  est_var <- matrix(NA,ncol=ntreat,nrow=2)
 } # end if

 # contrast matrix
 cont.mat <- matrix(c(diag(1,nrow=ntreat,ncol=ntreat),rep(-1,times=ntreat)),
                 nrow=ntreat,ncol=ntreat+1)


 # start simulations
 for (i in 1:nsim){

  # generate bivariate normal means
  tmean.early <- NULL
  tmean1.final <- NULL

  for (j in 1:ntreat){

   treatj <- simeans.binormal(n=nsamp[1],means=c(early[j],final[j]),vars=c(1,1),corr=corr)

   # interim
   tmean.early <- append(tmean.early,treatj$samp1)

   # final
   tmean1.final <- append(tmean1.final,treatj$samp2)

  } # end j

  # append control treatment
  ctmean1.final <- matrix(c(tmean1.final,sqrt(1/nsamp[1])*rnorm(1)),nrow=ntreat+1,ncol=1)
  ctmean.early <- matrix(c(tmean.early,sqrt(1/nsamp[1])*rnorm(1)),nrow=ntreat+1,ncol=1)

  # test statistic
  Etest.stat <- cont.mat%*%(ctmean.early*sqrt(nsamp[1]/2))
  Itest.stat <- cont.mat%*%(ctmean1.final*sqrt(nsamp[1]/2))
  Etest.stat <- c(as.numeric(Etest.stat))

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

    # sample size for stage 2
    if(real==TRUE){
     nsamp2 <- floor(nsamp_stage2/(sum(interim.select$select)+1))
    } else {
     nsamp2 <- nsamp[2]
    } # end if

    # final treatment means
    tmean2.final <- c(final,0)+sqrt(1/nsamp2)*rnorm(ntreat+1)
    Fz <- cont.mat%*%(tmean2.final*sqrt(nsamp2/2))
    Ftest.stat <- interim.select$z+Fz

    # second stage p-values
    dunnett2 <- dunnett.test(Ftest.stat,select=interim.select$select)

    # combination test
    weight <- nsamp[1]/(nsamp[1]+nsamp[2])
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


  # if estimate is true
  if (estim==TRUE){

   if (sum(interim.select$select)>0){

    # data storage
    sum1treat_1 <- sum1treat_1+nsamp[1]*tmean1.final*as.numeric(interim.select$select)
    sum1treat_2 <- sum1treat_2+nsamp2*tmean2.final[1:ntreat]*as.numeric(interim.select$select)
    sumN_1 <- sumN_1+nsamp[1]*as.numeric(interim.select$select)
    sumN_2 <- sumN_2+nsamp2*as.numeric(interim.select$select)
    sum2treat_1 <- sum2treat_1+nsamp[1]*((tmean1.final*as.numeric(interim.select$select))^2)
    sum2treat_2 <- sum2treat_2+nsamp2*((tmean2.final[1:ntreat]*as.numeric(interim.select$select))^2)
    # save results
    est_mean[1,] <- sum1treat_1/sumN_1
    est_mean[2,] <- sum1treat_2/sumN_2
    est_N[1,] <- sumN_1
    est_N[2,] <- sumN_2
    est_var[1,] <- (sum2treat_1-sumN_1*((est_mean[1,])^2))/(sumN_1-1)
    est_var[2,] <- (sum2treat_2-sumN_2*((est_mean[2,])^2))/(sumN_2-1)

   } # end if

  } # end if estim

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
 if (estim==TRUE){
  colnames(est_mean) <- as.character(1:ntreat)
  colnames(est_N) <- as.character(1:ntreat)
  colnames(est_var) <- as.character(1:ntreat)
  rownames(est_mean) <- c("Stage1","Stage2")
  rownames(est_N) <- c("Stage1","Stage2")
  rownames(est_var) <- c("Stage1","Stage2")
  list(count.total=count_total,select.total=select_total,
     reject.total=reject_total,sim.reject=sim.reject,
     est.mean=est_mean,est.N=est_N,est.var=est_var)
 } else {
  list(count.total=count_total,select.total=select_total,
     reject.total=reject_total,sim.reject=sim.reject)
 } # end if

} # end of function asd.sim

