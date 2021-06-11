require(MCMCpack)
require(sn)
require(combinat)
require(MASS)
### MAIN
fmm.overlap=function(y,K,maxT=5000,seed=1,method)
{
  ptm <- proc.time()
  set.seed(seed)
  d=ncol(y) ## number of events
  n=nrow(y) ## number of units
  K_star=2^K ## number of clusters
  u=z_ext(0:1,K) ## connection matrix
  prim=which(rowSums(u)==1) ### quali sono primary
  ## corrispondenza fra i primary e l'ordinamento k=1,2,3,...,K
  which.prim=rep(0,K) 
  for(k in 1:K) which.prim[k]=which(u[prim[k],]==1)
  ###
  z=matrix(0,n,K) ## original allocation vectors
  z_star=matrix(0,n,K_star) ## multiple allocation vectors
  
  
  # initialize alpha_star
  alpha_star.chain=matrix(0,maxT,K_star) ## canvas
  alpha_star.prior=rep(1,K_star)
  alpha_star.prior[prim]=K_star
  alpha_star=rdirichlet(1,alpha_star.prior)
  
  # initialize p
  
  pi.greco=matrix(0,K,d) ## matrix of attendance probs
  pi.greco.chain=array(0,c(maxT,K,d)) ## canvas
  pi.greco.prior=array(rep(1,K*d*2),c(K,d,2)) ## prior on attendance probs
  for (j in 1:d) for(k in 1:K) pi.greco[k,j]=rbeta(1,1,1) ## initialize with priors
  index=rep(0,n) 
  for (i in 1:n) {index[i]=sample(1:K_star,1,prob=alpha_star)
  z_star[i,index[i]]=1}
  z=u[index,]
  s=array(0,c(d,n,K)) ## auxiliary variable 
  for(j in 1:d){
    for (i in 1:n){
      if(sum(z[i,])<=1) s[j,i,]=z[i,] else s[j,i,which.max(z[i,]*pi.greco[,j])]=1
    }
  }
  
  prob.z_star.chain=array(0,c(maxT,n,K_star)) ### canvas for post.prob of allocation
  z_star.chain=array(0,c(maxT,n,K_star))
  loglik=rep(0,maxT)
  for (tt in 1:maxT) {
    
    ## Step 1. Sample z
    prob.z_star=matrix(0,n,K_star)
    
    for(h in 1:K_star){
      f.z.y=rep(0,n)
      for(j in 1:d) {
        prbl=ifelse(sum(u[h,])==0,0.00000001,psi.fun(pi.greco[,j],u[h,],method))
        f.z.y=f.z.y+apply(as.matrix(y[,j]),1,dbinom,1,prob=prbl,log=TRUE)
      }
      prob.z_star[,h]=log(alpha_star[h])+f.z.y
    }
    prob.z_star=exp(prob.z_star)
    prob.z_star=prob.z_star/rowSums(prob.z_star)
    if(any(is.na(prob.z_star))==1) print("FLAG: NA on prob.z")
    prob.z_star=ifelse(is.na(prob.z_star),1/K_star,prob.z_star)
    z_star=t(apply(prob.z_star,1,function(x) {rmultinom(1,1,prob=x)}))
    # z_star=true.param$z_star
    index=apply(z_star,1,which.max)
    
    
    ## Step 2. Allocate units according to sampled values
    z=u[index,]
    
    ### Step 6. compute s
    for(j in 1:d){
      for (i in 1:n){
        if(sum(z[i,])<=1){
          s[j,i,]=z[i,]
        }else{
          hlp=rep(0,K)
          if (method=="min") hlp[which.min(pi.greco[,j]^z[i,])]=1 else hlp[which.max(z[i,]*pi.greco[,j])]=1
          s[j,i,]=hlp
        }
      }
    }
    
    ## Step 3. Retrieve n_star
    n_star=colSums(z_star)
    
    # Step 4. Sample alpha_star
    alpha_star=rdirichlet(1,n_star+alpha_star.prior)
    # alpha_star=true.param$alpha_star
    
    # Step 5. Sample pi.greco
    for(k in 1:K){
      for (j in 1:d){
        pi.greco[k,j]<-rbeta(1,y[,j]%*%s[j,,k]+pi.greco.prior[k,j,1],sum(s[j,,k])-y[,j]%*%s[j,,k]+pi.greco.prior[k,j,2])
      }
    }
    # pi.greco=true.param$pi.greco
    
    # # ## Online Label Switching reassignment
    # criter=rowSums(pi.greco)
    # # # criter=alpha_star[prim]
    # relab=order(criter)
    # buff.u=u
    # # buff.pi=pi.greco
    # # # buff.alpha_star=alpha_star
    # for(k in 1:K){
    #   buff.u[,k]=u[,relab[k]]
    #   # buff.pi[k,]=pi.greco[relab[k],]
    #   # buff.alpha_star[prim[which.prim[k]]]=alpha_star[prim[relab[k]]]
    # }
    # u=buff.u
    # # # pi.greco=buff.pi
    # # # alpha_star=buff.alpha_star
    
    temp_sum=0
    for(j in 1:d){
      prbl=rep(0,n)
      for(i in 1:n){
        prbl[i]=ifelse(index[i]==1,0.00000001,psi.fun(pi.greco[,j],u[index[i],],"max"))
        temp_sum=temp_sum+y[i,j]*log(prbl[i])+(1-y[i,j])*log(1-prbl[i])
      }
    }

    
    ## saving chains into canvases
    prob.z_star.chain[tt,,]=prob.z_star
    pi.greco.chain[tt,,]=pi.greco
    alpha_star.chain[tt,]=alpha_star
    loglik[tt]=temp_sum
    ## print iteration number
    # print(tt)
  }
  
  out=(list(prob.z_star.chain=prob.z_star.chain,
            pi.greco.chain=pi.greco.chain,
            alpha_star.chain=alpha_star.chain,
            loglik=loglik,
            tempo=proc.time()-ptm))
  
}

########### Functions
## combining scheme
psi.fun<-function(x,u,method){
  if(method=="min") return(min(x^u)) else return(max(x*u))
}


### computes u connection matrix
z_ext <-function(x,nfac){
  ## x    : punti di quadratura o pesi
  ## nfac : numero di fattori
  nq <- length(x)                           # nr punti di quadratura
  zx <- hcube(rep(nq,nfac))                 # calcola tutte le possibili disposizioni con ripetizione di 8 elementi presi a gruppi di 3
  zx <- zx[,dim(zx)[2]:1]                   # zx contiene tutte le disposizioni "rovesciate"
  z2 <- matrix(x[zx],dim(zx)[1],dim(zx)[2]) # in corrispondenza di ciascuna posizione vado a selezionare il nodo corrispondente. In questo modo
  # ottengo tutte le possibili combinazioni tra quadrature per ogni fattore latente.
  return(z2)
}