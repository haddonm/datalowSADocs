

library(datalowSA)

# logpar=logparam; indat=hake; schaefer=TRUE; year="year";cpue="cpue"

simpspmL <- function(logpar,indat,schaefer=TRUE,year="year",cpue="cpue") { 
  ce <- indat[,cpue]
  nyrs <- length(indat[,year])
  biom <- numeric(nyrs+1)
  predCE <- numeric(nyrs)
  predC <- numeric(nyrs)
  r <- max(exp(logpar[1]),0.01)
  K <- exp(logpar[2])
  FF <- logpar[5:(5+nyrs)]
  biom[1] <- K
  if (logpar[3] != 0) biom[1] <- pars[2] * pars[3]
  #p is the location of mode parameter 1 = Schaefer, 1e-8 ~ Fox model
  if(schaefer)  p <- 1 else p <- 1e-8 
  for (index in 2:(nyrs+1)) {
    H <- 1/(1 + exp(-FF[index-1]))
    Bt <- biom[index-1]
    predC[index-1] <- H * Bt
    biom[index] <- max(Bt + ((r/p)*Bt*(1-(Bt/K)^p)-predC[index-1]),40)
  } 
  biom <- biom[1:nyrs]
  #for (i in 1:nce) {
    pick <- which(ce > 0)
    qcontrib <- log(ce[pick]/biom[pick]) # log(CE/Bt) each year's contrib to q
    q1 <- exp(mean(qcontrib,na.rm=TRUE))
    predCE[pick] <- biom[pick] * q1
 # }
  ans <- cbind(predC,rep(NA,nyrs))
  ans[pick,2] <- predCE
  colnames(ans) <- c("predC","predCE")
  return(ans)
}

hake <- read.table("C:/A_Mal/Rcodeold/TMB/tmb/LectC2.dat", header=TRUE)
names(hake) <- c("year", "catch", "cpue")

nyr <- length(hake$catch)
initFF <- rep(-2,nyr)

logparam <- c(-1.0,7.9,0,0.2,initFF)

simpspmL(logpar=logparam,indat=hake)





