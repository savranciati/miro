### Supreme Court voting: data from Doreian and Fujimoto (2003)
### Nine Justices:
# Covariates:
# year <- year they joined the Supreme Court
# topic <-  macro-category for the decisions taken by the SC
# Order of Justices in the data (rows) Breyer (Br), Ginsburg (Gi), Souter (So), Stevens (St), O'Connor (OC), Kennedy (Ke), Rehnquist (Re), Scalia (Sc) and Thomas (Th)
n <- 9
d <- 26
data <- matrix(0,n,d)
justices.names <- rownames(data) <- c("Breyer","Ginsburg",
                                      "Souter","Stevens",
                                      "O'Connor","Kennedy",
                                      "Rehnquist","Scalia","Thomas")
year <- c(1994,1993,
          1990,1975,
          1981,1988,
          1972,1982,1991)
topic <- c("Presidential Election",
           rep("Criminal law",5),
           rep("Federal authority",6),
           rep("Civil rights",3),
           rep("Immigration law",4),
           rep("Speech and Press",5),
           rep("Labor and Properties",2))
### "-" coded with 0, voted in the minority
### "+" coded with 1, voted in the majority
### "0" coded with 1, for simplicity
data[1,] <- c(0,
              1,1,1,0,1,
              0,1,0,1,0,0,
              1,0,1,
              1,1,1,0,
              1,1,1,1,0,
              0,0)
data[2,] <- c(0,
              1,1,1,0,1,
              0,1,0,1,0,1,
              1,0,1,
              1,1,1,0,
              1,1,0,1,0,
              0,0)  
data[3,] <- c(0,
              1,1,1,1,1,
              0,1,0,1,1,1,
              1,0,1,
              1,1,1,0,
              1,1,0,1,0,
              0,0)
data[4,] <- c(0,
              1,1,0,0,1,
              0,1,0,1,1,0,
              1,0,1,
              1,1,1,1,
              1,1,0,1,0,
              0,0)
data[5,] <- c(1,
              1,1,0,0,1,
              1,1,1,1,0,1,
              1,1,1,
              0,1,0,0,
              0,1,1,1,1,
              1,1)
data[6,] <- c(1,
              1,1,0,1,1,
              1,1,1,1,1,1,
              0,1,1,
              1,0,1,1,
              1,1,1,0,1,
              1,1)
data[7,] <- c(1,
              0,0,0,1,0,
              1,1,1,1,1,1,
              0,1,1,
              0,0,0,1,
              0,0,1,0,1,
              1,1)
data[8,] <- c(1,
             0,0,1,1,0,
             1,1,1,1,1,1,
             0,1,0,
             0,0,0,1,
             0,0,1,0,1,
             1,1)
data[9,] <- c(1,
              0,0,1,1,0,
              1,1,1,1,1,1,
              0,1,0,
              0,0,0,1,
              0,0,1,0,1,
              1,1)

source('manet.R')
source('mixt_Bern.R')
source('mixt_probit.R')
source('miro_event.R')

###############################################################
######################### mixt_probit ###############################
###############################################################

fake.X=matrix(0,n,2)
maxT=10000

require(dummies)
matW=dummy(topic,sep=".")
matW=cbind(rep(1,d),matW[,-1])

### mixtprobit
mod_mixtprobit <- mixt_probit(data,X=fake.X,W=matW,2,maxT=maxT)
burn.in=round(maxT*0.9,0)
post.prob.z=apply(mod_mixtprobit$prob.z.chain[burn.in:maxT,,],c(2,3),mean)
post.prob.z=post.prob.z/rowSums(post.prob.z)
clasf.mixtprobit=apply(post.prob.z,1,which.max)
table(clasf.mixtprobit)

n.par.mixtprobit=2+7*2
AICM.mixtprobit=-2*(mean(mod_mixtprobit$log.lik[burn.in:maxT],na.rm=TRUE)-var(mod_mixtprobit$log.lik[burn.in:maxT],na.rm=TRUE))
BIC_mcmc.mixtprobit=-2*max(mod_mixtprobit$log.lik[burn.in:maxT],na.rm=TRUE)+n.par.mixtprobit*log(n*d)
### the lower the better
AICM.mixtprobit
BIC_mcmc.mixtprobit


###############################################################
######################### mixtBern ###############################
###############################################################

maxT=10000

### mixtbern
mod_mixtbern <- fmm.no_over(data,3,maxT=maxT)
burn.in=round(maxT*0.9,0)
post.prob.z=apply(mod_mixtbern$prob.z.chain[burn.in:maxT,,],c(2,3),mean)
post.prob.z=post.prob.z/rowSums(post.prob.z)
clasf.mixtbern=apply(post.prob.z,1,which.max)
table(clasf.mixtbern)

n.par.mixtbern=3+26*3
AICM.mixtbern=-2*(mean(mod_mixtbern$loglik[burn.in:maxT],na.rm=TRUE)-var(mod_mixtbern$loglik[burn.in:maxT],na.rm=TRUE))
BIC_mcmc.mixtbern=-2*max(mod_mixtbern$loglik[burn.in:maxT],na.rm=TRUE)+n.par.mixtbern*log(n*d)
### the lower the better
AICM.mixtbern
BIC_mcmc.mixtbern


###############################################################
######################### manet ###############################
###############################################################

maxT=10000
### manet
mod_manet <- fmm.overlap(data,2,maxT=maxT,method="max")
# post
burn.in=round(maxT*0.9,0)
post.prob.z_star=apply(mod_manet$prob.z_star.chain[burn.in:maxT,,],c(2,3),mean)
post.prob.z_star=post.prob.z_star/rowSums(post.prob.z_star)
clasf.manet=apply(post.prob.z_star,1,which.max)
table(clasf.manet)

n.par.manet=4+26*2
AICM.manet=-2*(mean(mod_manet$loglik[burn.in:maxT],na.rm=TRUE)-var(mod_manet$loglik[burn.in:maxT],na.rm=TRUE))
BIC_mcmc.manet=-2*max(mod_manet$loglik[burn.in:maxT],na.rm=TRUE)+n.par.manet*log(n*d)
### the lower the better
AICM.manet
BIC_mcmc.manet


t(round(apply(mod_manet$pi.greco.chain[burn.in:maxT,,],c(2,3),mean),3))
t(round(apply(mod_manet$pi.greco.chain[burn.in:maxT,,],c(2,3),median),3))
t(round(apply(mod_manet$pi.greco.chain[burn.in:maxT,,],c(2,3),sd),3))

###############################################################
##################### miro_event ##############################
###############################################################
maxT=10000
### only event covariate: topic
matW=dummy(topic,sep=".")
matW=cbind(rep(1,d),matW[,-1])
mod_miro <- miro_mean(data,W=matW,K=2,maxT=maxT)

# post
burn.in=round(maxT*0.9,0)
post.prob.z_star=apply(mod_miro$prob.z_star.chain[burn.in:maxT,,],c(2,3),mean)
post.prob.z_star=post.prob.z_star/rowSums(post.prob.z_star)
clasf.miro=apply(post.prob.z_star,1,which.max)
table(clasf.miro)

n.par.miro <- 4+14
AICM.miro=-2*(mean(mod_miro$log.lik[burn.in:maxT],na.rm=TRUE)-var(mod_miro$log.lik[burn.in:maxT],na.rm=TRUE))
BIC_mcmc.miro=-2*max(mod_miro$log.lik[burn.in:maxT],na.rm=TRUE)+n.par.miro*log(n*d)
### the lower the better
AICM.miro
BIC_mcmc.miro
