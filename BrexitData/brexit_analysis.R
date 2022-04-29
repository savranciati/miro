library(dummies)
library(ggplot2)
library(reshape2)
library(forcats)
library(dplyr)
source('miro_wrapper.R')
bin_ayeno <- function(x){
  ### partial analysis
  # return(ifelse(x=="Aye_Vote", 1,0))
  ### full analysis
  return(ifelse(x==("Aye_Vote" | "No_Vote"), 1,0))
}

### Load/prepare data
load("Brexit.RData")
parties <- Brexit$Party
MPnames <- rownames(Brexit)
Leave <- Brexit$Leave
EE <- Brexit$ExpEntropy
Age <- Brexit$Age
Brex.data <- Brexit[,c(6:13)]
colnames(Brex.data) <- c("(B)_NoDeal","(D)_CommonMarket",
                         "(E)_EFTA_EEA","(J)_CustomsUnion",
                         "(K)_LaboursPlan","(L)_RevocAvoidNoDeal",
                         "(M)_ConfirmPublicVote","(O)_ManagedNoDeal")


{

  NoAbs_subset <- which((apply(Brex.data,1,function(x) sum(x=="absent"))==0)==T)
  apply(Brex.data[NoAbs_subset,],2,table)
  
  parties <- parties[NoAbs_subset]
  parties <- factor(parties)
  MPnames <- MPnames[NoAbs_subset]
  MPnames <- factor(MPnames)
  Leave <- Leave[NoAbs_subset]
  EE <- EE[NoAbs_subset]
  Age <- Age[NoAbs_subset]
  dt.Brexit <- Brex.data[NoAbs_subset,]
  # rm(Brex.data,Brexit, NoAbs_subset)
  
  ### Convert to binary
  # Original factors: 0 absent 1 Aye 2 No
  for(ii in 1:8){
    # 0 <- Aye , 1 <- No
    dt.Brexit[,ii] <- as.numeric(dt.Brexit[,ii])-2
    # 0 <- No , 1 <- Aye
    dt.Brexit[,ii] <- as.numeric(!dt.Brexit[,ii])
  }
  str(dt.Brexit)
  
  
}


### Create covariate for divisions (variable-specific)
divis.cov <- factor(c("ProBrexit",rep("AgainstBrexit",6),"ProBrexit"))
matW=dummy(divis.cov,sep=".")
matW=as.matrix(matW[,-1])
colnames(matW) <- c("ProBrexit")
## different coding for ProBrexit
matW[2:7] <- -1

### Create covariate for MPs (Leave)
matX <- cbind(Leave,EE)

### Fit model
n <- dim(dt.Brexit)[1]
d <- dim(dt.Brexit)[2]
K <- 2
maxT <- 10000
time.start <- Sys.time()
### with interaction
mod.brexit <- miro_int(as.matrix(dt.Brexit),X=matX,W=matW,K,prior.var=1,
                       maxT=maxT,seed=42, nullClust=F)

time.end <- Sys.time()

### post analysis
time.end-time.start
burn.in=round(maxT*0.5,0)
### Plot lik
plot(mod.brexit$log.lik[burn.in:maxT], type="l")

post.prob.z_star=apply(mod.brexit$prob.z_star.chain[burn.in:maxT,,],c(2,3),mean)
post.prob.z_star=post.prob.z_star/rowSums(post.prob.z_star)
clasf.miro=apply(post.prob.z_star,1,which.max)

table(clasf.miro)
table(clasf.miro, parties)

### relabel to match vectors
buff.clasf <- clasf.miro
clasf.miro[clasf.miro==1] <- 2
clasf.miro[buff.clasf==2] <- 1

table(clasf.miro, parties)

n.par.miro <- dim(mod.brexit$alpha_star.chain)[2]+K*(dim(matX)[2]+dim(matW)[2]+dim(matX)[2]*dim(matW)[2])
AICM.manetcov=-2*(mean(mod.brexit$log.lik[burn.in:maxT],na.rm=TRUE)-var(mod.brexit$log.lik[burn.in:maxT],na.rm=TRUE))
BIC_mcmc.manetcov=-2*max(mod.brexit$log.lik[burn.in:maxT],na.rm=TRUE)+n.par.miro*log(n*d)
### the lower the better
AICM.manetcov
BIC_mcmc.manetcov

