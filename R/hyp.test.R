hyp.test <-
function(comb.test,level=level,full.hyp=FALSE){

 # test level
 crit <- qnorm(1-level)

 # number of treatments is length of list
 treats <- length(comb.test$hyp.comb)

 # combined zscores
 zscores <- comb.test$zscores
 hyp.reject <- matrix(0,nrow=1,ncol=treats)
 colnames(hyp.reject) <- colnames(comb.test$zscores[[1]])
 rownames(hyp.reject) <- 1

 # rejection indicator
 hyp.reject <- matrix(0,nrow=1,ncol=treats)
 colnames(hyp.reject) <- colnames(zscores[[1]])
 rownames(hyp.reject) <- 1

 # unlist zscores
 vhyp.comb <- unlist(comb.test$hyp.comb)
 hyp.len <- rep(1:treats,times=choose(treats,1:treats))
 vzscores <- rep(unlist(comb.test$zscores),hyp.len)

 # matrix of zscores
 matzscores <- matrix(NA,nrow=length(vzscores)/treats,ncol=treats)
 matrejects <- matzscores
 for (i in 1:treats){
  matzscores[,i] <- vzscores[vhyp.comb==i]
  if (full.hyp==TRUE) {matrejects[,i] <- as.numeric(matzscores[,i]>=crit)}
  if (sum(as.numeric(matzscores[,i]>=crit))==length(matzscores[,i])){
   hyp.reject[i] <- 1
  } # end if
 } # end i
 colnames(matzscores) <- colnames(comb.test$zscores[[1]])
 colnames(matrejects) <- colnames(comb.test$zscores[[1]])
 rownames(matzscores) <- 1:(length(vzscores)/treats)
 rownames(matrejects) <- 1:(length(vzscores)/treats)
 
 # matrix of intersection hypotheses
 if (full.hyp==TRUE){
  hyp.names <- NULL
  for (i in 1:treats){
   hyp.names <- append(hyp.names, colnames(comb.test$hyp.comb[[i]]))
  } # end i
  vhypcomb <- rep(hyp.names,hyp.len)
  mathypcomb <- matrix(NA,nrow=length(vzscores)/treats,ncol=treats)
  for (i in 1:treats){
   mathypcomb[,i] <- vhypcomb[vhyp.comb==i]
  } # end i
  colnames(mathypcomb) <- colnames(comb.test$zscores[[1]])
  rownames(mathypcomb) <- 1:(length(vzscores)/treats)
 } # end full.hyp=TRUE

 # results of hypotheses tests
 if (full.hyp==TRUE){
  list(reject=hyp.reject,all.rejects=matrejects,all.hyp=mathypcomb)
 } else {
  list(reject=hyp.reject)
 } # end if

} # end of function hyp.test

