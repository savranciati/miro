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
### manet
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
##################### manet+cov ###############################
###############################################################
maxT=10000
### solo event covariate: topic
require(dummies)
source('manet+cov_mean_event.R')
matW=dummy(topic,sep=".")
matW=cbind(rep(1,d),matW[,-1])
mod_manetcov <- manet.cov_mean(data,W=matW,K=2,maxT=maxT)

# post
burn.in=round(maxT*0.9,0)
post.prob.z_star=apply(mod_manetcov$prob.z_star.chain[burn.in:maxT,,],c(2,3),mean)
post.prob.z_star=post.prob.z_star/rowSums(post.prob.z_star)
clasf.manetcov=apply(post.prob.z_star,1,which.max)
table(clasf.manetcov)

n.par.manetcov <- 4+14
AICM.manetcov=-2*(mean(mod_manetcov$log.lik[burn.in:maxT],na.rm=TRUE)-var(mod_manetcov$log.lik[burn.in:maxT],na.rm=TRUE))
BIC_mcmc.manetcov=-2*max(mod_manetcov$log.lik[burn.in:maxT],na.rm=TRUE)+n.par.manetcov*log(n*d)
### the lower the better
AICM.manetcov
BIC_mcmc.manetcov

t(round(apply(mod_manetcov$gamma.chain[,burn.in:maxT],1,mean),3))
t(round(apply(mod_manetcov$gamma.chain[,burn.in:maxT],1,median),3))
t(round(apply(mod_manetcov$gamma.chain[,burn.in:maxT],1,sd),3))

par(mfrow=c(3,3))
plot(mod_manetcov$gamma.chain[1,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[1,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[2,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[2,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[3,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[3,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[4,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[4,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[5,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[5,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[6,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[6,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[7,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[7,burn.in:maxT]),col=2,lwd=2)
par(mfrow=c(1,1))

par(mfrow=c(3,3))
plot(mod_manetcov$gamma.chain[8,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[8,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[9,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[9,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[10,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[10,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[11,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[11,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[12,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[12,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[13,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[13,burn.in:maxT]),col=2,lwd=2)
plot(mod_manetcov$gamma.chain[14,burn.in:maxT],type="l")
abline(h=mean(mod_manetcov$gamma.chain[14,burn.in:maxT]),col=2,lwd=2)
par(mfrow=c(1,1))


plot(seq(1,7),t(round(apply(mod_manetcov$gamma.chain[,burn.in:maxT],1,mean),3))[1:7],col=1, type="b",ylim=c(-5,5))
lines(seq(1,7),t(round(apply(mod_manetcov$gamma.chain[,burn.in:maxT],1,mean),3))[8:14],col=2, type="b")
lines(seq(1,7),t(round(apply(mod_manetcov$gamma.chain[,burn.in:maxT],1,mean),3))[1:7]*0.5+t(round(apply(mod_manetcov$gamma.chain[,burn.in:maxT],1,mean),3))[8:14]*0.5,col=4, type="b")

differences <- abs(t(round(apply(mod_manetcov$gamma.chain[,burn.in:maxT],1,mean),3))[1:7]-t(round(apply(mod_manetcov$gamma.chain[,burn.in:maxT],1,mean),3))[8:14])
diff.data <- data.frame(differences,topic=sort(unique(topic)))
buff <- data.frame(topic)
topic2 <- merge(buff,diff.data,by="topic")
sorted.effects <- topic2[order(topic2$differences),]

#### PLOT manet+cov results
load("~/Desktop/temp manet+cov/dati reali manet+cov/Supreme Court dataset/supreme.RData")
require(ggplot2)
require(reshape2)
library(forcats)
detailed.topic <- c("PresElect.",
                    "IllegalSrc.1",
                    "IllegalSrc.2",
                    "IllegalSrc.3",
                    "SeatBelts",
                    "StayExec",
                    "Federalism",
                    "CleanAirAct",
                    "CleanWater",
                    "CannabisHealt",
                    "UnitedFood",
                    "NYTcopyrights",
                    "VotingRights",
                    "T.VI.Disabil",
                    "PGAvsHandicap",
                    "Immgr.Jurisdic",
                    "DeportCrimeAli",
                    "DetainCrimeAli",
                    "Citizenship",
                    "LegalAidPoor",
                    "Privacy",
                    "FreeSpeech",
                    "CampaignFinance",
                    "TobaccoAds",
                    "LaborRights",
                    "PropertyRights")
topic.merge <- data.frame(topic,detailed.topic)
topic3 <- unique(merge(topic.merge,sorted.effects,by="topic",all=FALSE))
rownames(topic3) <- NULL
colnames(data) <- factor(topic)
melted.data <- melt(data)
names(melted.data) <- c("Justice","topic","Voting")
levels(melted.data$topic)
levels(topic3$topic)
topic3$ordine <- c(rep(4,3),
                   rep(2,5),
                   rep(3,6),
                   rep(5,4),
                   rep(7,2),
                   rep(1,1),
                   rep(6,5))
topic3 <- topic3[order(topic3$ordine),]
rep.topic3 <- topic3[rep(seq_len(nrow(topic3)), each=9),]
rownames(rep.topic3) <- NULL
aug.data <- data.frame(melted.data,rep.topic3)
aug.data$topic.1 <- NULL
aug.data$ordine <- NULL
aug.data$det <- factor(aug.data$detailed.topic,ordered = TRUE, levels = levels(aug.data$detailed.topic))
aug.data$just <- factor(aug.data$Justice,ordered = TRUE, levels = levels(aug.data$Justice))
aug.data$vot <- factor(aug.data$Voting)

main_plot=ggplot(aug.data, aes(x=reorder(det,differences),
                     y=reorder(just,differences),
                     fill=reorder(vot,differences)))+
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, size=12),
        axis.text.y = element_text(hjust = 1, size=11))+
  scale_fill_manual(values=c("black","white"))+
  geom_tile()+xlab("Topic") + ylab("Justice")+
  labs(fill = "Voting")+
  geom_hline(yintercept=4.5, lwd=1.25, linetype=1, color = "white")+
  geom_hline(yintercept=6.5, lwd=1.25, linetype=1, color = "white")
main_plot

probs1 <- (t(round(apply(mod_manetcov$gamma.chain[,burn.in:maxT],1,mean),3)))[1:7]
probs2 <- (t(round(apply(mod_manetcov$gamma.chain[,burn.in:maxT],1,mean),3)))[8:14]
propro <- c(probs1,probs2)
cluster <- factor(c(rep(1,7),rep(2,7)))

a <- cumsum(c(3,5,5,4,6,1,2))
midpoints <- c(1.5,a[-length(a)] + diff(a)/2)
midpoints <- (midpoints/26)*7
midpoints <- round(midpoints,3)
probs <- data.frame(topic=rep(diff.data$topic,2),
                    propro,
                    ordine=rep(c(1,3,5,4,7,6,2),2),
                    cluster=cluster,
                    ticks=midpoints)
                    # ticks=round((cumsum(c(3,5,5,4,6,1,2))/26)*7,2))
probs$sorted.prob <- probs$propro[order(probs$ordine)]
probs$sorted.prob[1:7] <- probs$propro[order(probs$ordine[1:7])]
probs$sorted.prob[8:14] <- probs$propro[order(probs$ordine[8:14])+7]
etic <- as.character(unique(probs$topic[order(probs$ordine)]))

ggplot(probs, aes(x=ticks,y=sorted.prob, group=cluster)) +
         geom_line(aes(linetype = cluster))+
         geom_point(aes(shape=cluster))+
      theme(axis.title.y = element_text(size=15,margin=margin(1,25,1,1)),
        axis.text.x = element_text(angle = 90, hjust = 1, size=13),
      axis.text.y = element_text(hjust = 1, size=13))+
      xlab("Category") + ylab(expression(gamma))+
      scale_x_discrete(limits=probs$ticks[1:7],
                       labels=etic,
                       expand = expand_scale(mult=c(0.09,0.05)))+
      scale_y_continuous(limits=c(-4.2,4.2),
                         breaks=seq(-4.2,4.2,length.out=7),
                         expand=c(0,0.03))
      



      
  
df <- as.data.frame(data)
library(igraph)
temp=graph.incidence(t(df),mode=c("all"))

vert.shp <- c(rep("square",26),rep("circle",9))
vert.col <- rep(0,35)
vert.frm <- rep(0,35)
vert.size <- rep(0,35)
vert.labcol <-  rep(rgb(0,0,0,1),35)
vert.lab <- c(rep("",26),
              substr(justices.names,1,3))

vert.labfont <- c(rep(1,26),rep(2,9))
vert.col[1:26] <- rep(rgb(0.5,0.5,0.5,0.5),26)
vert.size[1:26] <- 9
vert.size[27:35] <- 16
vert.labcol[1:35]<-rep(rgb(0,0,0,1),35)
vert.frm[1:35] <- rep(rgb(0,0,0,1),35)

set.seed(100)
png("judg.png",height=1000,width=1300,units="px",res=96)
plot.igraph(temp,layout=layout.auto,
            vertex.color=vert.col,
            vertex.shape=vert.shp,
            vertex.frame.color=vert.frm,
            vertex.size=vert.size,
            vertex.label=vert.lab,
            vertex.label.color=vert.labcol,
            label.font=vert.labfont,
            vertex.label.cex=1.2)
dev.off()



















