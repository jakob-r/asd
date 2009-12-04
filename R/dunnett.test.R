dunnett.test <-
function(Z=Z,select=rep(1,length(Z))){

 # number of treatments is length of Z vector
 treats <- length(Z)

 # intersection hypotheses indicator
 hyp.comb <- list(NULL); nhyp.comb <- vector(length=treats)
 for (i in 1:treats){
  comb.dat <- combn(1:treats,i)
  hyp.comb[[i]] <- comb.dat
  nhyp.comb[i] <- ncol(hyp.comb[[i]])
  rownames(hyp.comb[[i]]) <- 1:i
  hypcol <- NULL
  for (j in 1:nhyp.comb[i]){
   ihypcol <- paste("H",paste(hyp.comb[[i]][,j],collapse=""),sep="")
   hypcol <- append(hypcol,ihypcol)
  } # end j
  colnames(hyp.comb[[i]]) <- hypcol
 } # end i

 # test results
 pdunnett.test <- list(NULL)
 zscores <- list(NULL)

 # dunnett function
 int_dunnett <- function(x,z,k){
  ((pnorm(sqrt(2)*z+x))^k)*dnorm(x)
 }

 # evaluate integral 
 for (i in 1:treats){

  # i = number of elements in intersection hypothesis
  ptest <- NULL

  # set Z to -Inf if not selected
  if(select[i]==0){
   Z[i] <- -Inf
  } # end if

  for (j in 1:nhyp.comb[i]){

   # calculate correct k
   kselect <- sum(select[c(hyp.comb[[i]][,j])])

   # select maximum 
   Zmax <- max(Z[c(hyp.comb[[i]][,j])])

   # calculate integral for z=Zmax and k=kselect
   if (kselect==0){
    dunnet_integral <- 0
    F_Zmax <- 1
   } else {
    dunnett_integral <- integrate(int_dunnett,-Inf,Inf,z=Zmax,k=kselect)
    F_Zmax <- 1-dunnett_integral$value
   } # end if

   ptest <- append(ptest,F_Zmax)
   
  } # end j

  # store results
  pdunnett.test[[i]] <- matrix(ptest,nrow=1,ncol=length(ptest))
  colnames(pdunnett.test[[i]]) <- colnames(hyp.comb[[i]])
  rownames(pdunnett.test[[i]]) <- 1
  zscores[[i]] <- qnorm(1-pdunnett.test[[i]])

 } # end i

 # results of Dunnett's tests
 list(pvalues=pdunnett.test,zscores=zscores,hyp.comb=hyp.comb)

} # end of function dunnett.test

