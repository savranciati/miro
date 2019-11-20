require(combinat)
require(sn)
require(MCMCpack)
require(mvtnorm)
require(truncnorm)
### Model n.3 - Overlapping cluster
### cluster-dependent coeff. for actors and events covariates

## Provide:
# y, data matrix "n x d" with n actors and d events
# W, design matrix "n x Q" where Q number of event covariates +1 for the intercept
# to be passed as dummy variables (1st column is intercept)
# K, number of overlapping clusters
# prior.var, variance for the prior on the regression coefficients \gamma and \gamma (default is 10)
# maxT, number of MCMC iterations (default is 5000)
# seed
# method, "min" or "max", default is "min"
manet.cov_mean<-function(y,W,K,prior.var=5,maxT=5000,seed=42)
{
  set.seed(seed)
  d=ncol(y) ## number of events
  n=nrow(y) ## number of units
  Q=ncol(W) ## number of actor covariates (intercept included)
  NN=n*d ## stacked vectors length
  K_star=2^K ## number of clusters
  u=z_ext(0:1,K) ## connection matrix
  kron.u=u/rowSums(u) ## connection matrix for kronecker product
  
  ## stacked matrices and vectors
  stack.y=c(y) # stacking by events
  stack.W=W
  rows=rep(1:nrow(stack.W),cbind(1:nrow(stack.W),n)[,2])
  stack.W=stack.W[rows,]
  
  ### Initialize parameters mixture weights
  # alpha_star.chain=matrix(0,maxT,K_star) ## canvas
  # alpha_star=rep(0,K_star) ## star vector 
  # alpha_star.prior=rep(1,K_star)
  # alpha_star=rdirichlet(1,alpha_star.prior)
  prim=which(rowSums(u)==1) ### quali sono primary
  alpha_star.chain=matrix(0,maxT,K_star) ## canvas
  alpha_star.prior=rep(1,K_star)
  alpha_star.prior[prim]=K_star
  alpha_star=rdirichlet(1,alpha_star.prior)
  
  ## Initialize allocation variables z
  z=matrix(0,n,K) ## original allocation vectors
  z_star=matrix(0,n,K_star) ## multiple allocation vectors
  index=rep(0,n)
  for (i in 1:n){
    index[i]=sample(1:K_star,1,prob=alpha_star)
    z_star[i,index[i]]=1
  }
  z=u[index,]
  
  ## Build design matrix
  ## Compute linear predictor and probs
  ## initialize regr.coeff from prior
  design=array(0,c(n*d,Q*K,K_star))
  lin.pred=matrix(0,n*d,K_star)
  prior.m.coeff=rep(0,Q*K) 
  prior.var.coeff=prior.var*diag(Q*K)
  gamma.chain=matrix(0,Q*K,maxT)
  temp.coeff=c(rmvnorm(1,mean=prior.m.coeff,sigma=prior.var.coeff))
  for(h in 1:K_star){
    design[,,h]=kronecker(t(kron.u[h,]),stack.W)
    lin.pred[,h]=design[,,h]%*%temp.coeff
  }
  probs=pnorm(lin.pred)
  probs[is.na(probs)]=0.00001
  probs[probs==0]=0.00001
  probs[probs==1]=0.99999
  
  ### canvas for post.prob of allocation
  prob.z_star.chain=array(0,c(maxT,n,K_star))
  z_star.chain=array(0,c(maxT,n,K_star))
  
  ### canvas for log.likelihood
  flag.na.z=log.lik=rep(0,maxT)
  
  ##### Start MCMC
  for (tt in 1:maxT){
    
    ## Step 1. Sample z
    prob.z_star=matrix(0,n,K_star)
    for(h in 1:K_star){
      f.y.z=stack.y*log(probs[,h])+(1-stack.y)*log(1-probs[,h])
      prob.z_star[,h]=log(alpha_star[h])+rowSums(matrix(f.y.z,n,d,byrow=FALSE))
    }
    prob.z_star=exp(prob.z_star)
    prob.z_star=prob.z_star/rowSums(prob.z_star)
    if(any(is.na(prob.z_star))==1) flag.na.z[tt]=1
    prob.z_star=ifelse(is.na(prob.z_star),1/K_star,prob.z_star)
    prob.z_star.chain[tt,,]=prob.z_star
    z_star=t(apply(prob.z_star,1,function(x) {rmultinom(1,1,prob=x)}))
    index=apply(z_star,1,which.max)
    n_star=colSums(z_star)
    # index=true.params$index
    # n_star=table(index)

    full.index=cbind(seq_along(rep(index,d)),rep(index,d))
    
    ## Step 2. Retrieve n_star and Sample alpha_star
    alpha_star=rdirichlet(1,n_star+alpha_star.prior)
    alpha_star.chain[tt,]=alpha_star
  
    design.z=matrix(0,NN,Q*K)
    ## Step 3. build design matrix conditional on allocation
    for(pp in 1:(Q*K)){
      buff=design[,pp,]
      design.z[,pp]=buff[full.index]
    }
    
    if(n_star[1]>0) rem_zero=(full.index[,2]==1) else rem_zero=rep(FALSE,NN)
  
    ## Step 4. Sample regr.coeff; compute lin.pred and probs
    out.coeff=probit_sampler(y=stack.y[!rem_zero],
                            cov.mat=design.z[!rem_zero,],
                            buff.coeff=temp.coeff,
                            prior.m.coeff=prior.m.coeff,
                            prior.var.coeff=prior.var.coeff)
    temp.coeff=out.coeff
    gamma.chain[,tt]=temp.coeff
    
    ## Step 5. update preds and pnorms
    for(h in 1:K_star){
      lin.pred[,h]=design[,,h]%*%temp.coeff
    }
    probs=pnorm(lin.pred)
    probs[is.na(probs)]=0.00001
    probs[probs==0]=0.00001
    probs[probs==1]=0.99999
    
    ### compute log.lik
    log.lik[tt]=sum(stack.y*log(probs[full.index])+(1-stack.y)*log(1-probs[full.index]))
    
    ## print iteration number
    print(tt)
  }
  
  out=(list(prob.z_star.chain=prob.z_star.chain,
            alpha_star.chain=alpha_star.chain,
            gamma.chain=gamma.chain,
            last.cl=index,
            log.lik=log.lik,
            flag.na.z=flag.na.z))
}

########### Functions
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
#### basic sampler inside tt loop from iterative sampling (Holmes, 2006)
probit_sampler<-function(y,cov.mat,buff.coeff,prior.m.coeff,prior.var.coeff){
  B=0
  aux=rep(0,nrow(cov.mat))
  post.var=solve(t(cov.mat)%*%cov.mat+solve(prior.var.coeff))
  lower.post.var=t(chol(post.var))
  S=post.var%*%t(cov.mat)
  for(i in 1:nrow(cov.mat)){
    m=cov.mat[i,]%*%buff.coeff
    if(y[i]==1) aux[i]=rtruncnorm(1,a=0,b=Inf,mean=m,sd=1) else aux[i]=rtruncnorm(1,a=-Inf,b=0,mean=m,sd=1)
    B=B+aux[i]*S[,i]
  }
  temp=rmvnorm(1,prior.m.coeff,prior.var.coeff)
  return(temp.coeff=B+lower.post.var%*%t(temp))
}
