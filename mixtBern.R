require(MCMCpack)
require(sn)
require(MASS)
### MAIN
fmm.no_over=function(y,K,maxT=5000,seed=1)
{
  ptm <- proc.time()
  set.seed(seed)
  d=ncol(y) ## number of events
  n=nrow(y) ## number of units
  z=matrix(0,n,K) ## original allocation vectors
  
  
  # initialize alpha
  alpha.chain=matrix(0,maxT,K) ## canvas
  alpha.prior=rep(1,K) ## prior 
  alpha=rep(1/K,K) ## vector
  
  # initialize p
  
  pi.greco=matrix(0,K,d) ## matrix of attendance probs
  pi.greco.chain=array(0,c(maxT,K,d)) ## canvas
  pi.greco.prior=array(rep(1,K*d*2),c(K,d,2)) ## prior on attendance probs
  for (j in 1:d) for(k in 1:K) pi.greco[k,j]=rbeta(1,1,1) ## initialize with priors
  index=rep(0,n) 
  for (i in 1:n) {index[i]=sample(1:K,1,prob=alpha)
  z[i,index[i]]=1}
  
  prob.z.chain=array(0,c(maxT,n,K)) ### canvas for post.prob of allocation
  z.chain=array(0,c(maxT,n,K))
  loglik=rep(0,maxT)
  for (tt in 1:maxT) {
    
    ## Step 1. Sample z
    prob.z=matrix(0,n,K)
    
    for(k in 1:K){
      f.z.y=rep(0,n)
      for(j in 1:d) {
        f.z.y=f.z.y+apply(as.matrix(y[,j]),1,dbinom,1,prob=pi.greco[k,j],log=TRUE)
      }
      prob.z[,k]=log(alpha[k])+f.z.y
    }
    prob.z=exp(prob.z)
    prob.z=prob.z/rowSums(prob.z)
    prob.z=ifelse(is.na(prob.z),1/K,prob.z)
    z=t(apply(prob.z,1,function(x) {rmultinom(1,1,prob=x)}))
    
    ## Step 3. Retrieve z and n_k
    n_k=colSums(z)
    
    # Step 5. Sample alpha_star
    alpha=rdirichlet(1,n_k+alpha.prior)
    
    
    # Step 6. Sample pi.greco
    temp_sum=0
    for(k in 1:K){
      for (j in 1:d){
        pi.greco[k,j]<-rbeta(1,y[,j]%*%z[,k]+pi.greco.prior[k,j,1],sum(z[,k])-y[,j]%*%z[,k]+pi.greco.prior[k,j,2])
        temp_sum=temp_sum+sum(y[z[,k]==1,j]*log(pi.greco[k,j]))+sum((1-y[z[,k]==1,j])*log(1-pi.greco[k,j]))
      }
    }
    

    ## saving chains into canvases
    prob.z.chain[tt,,]=prob.z
    pi.greco.chain[tt,,]=pi.greco
    alpha.chain[tt,]=alpha
    z.chain[tt,,]=z
    loglik[tt]=temp_sum
    ## print iteration number
    # print(tt)
  }
  
  out=(list(prob.z.chain=prob.z.chain,
            z.chain=z.chain,
            pi.greco.chain=pi.greco.chain,
            alpha.chain=alpha.chain,
            loglik=loglik,
            tempo=proc.time()-ptm))
  
}