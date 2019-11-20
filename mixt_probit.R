require(combinat)
require(sn)
require(MCMCpack)
require(mvtnorm)
require(truncnorm)
### Mixture of Probit Regression Models
### cluster-dependent coeff. for actors and events covariates

## Provide:
# y, data matrix "n x d" with n actors and d events
# X, design matrix "n x L" where L number of actor covariates +1 for the intercept
# W, design matrix "d x Q" where Q number of event covariates
# to be passed as dummy variables (1st column is intercept)
# K, number of clusters
# prior.var, variance for the prior on the regression coefficients \gamma and \gamma (default is 5)
# maxT, number of MCMC iterations (default is 5000)
# seed
mixt_probit<-function(y,X,W,K,prior.var=2,maxT=5000,seed=42)
{
  set.seed(seed)
  d=ncol(y) ## number of events
  n=nrow(y) ## number of units
  L=ncol(X) ## number of actor covariates (intercept included)
  Q=ncol(W) ## number of actor covariates
  NN=n*d ## stacked vectors length
  
  ## stacked matrices and vectors
  stack.y=c(y) # stacking by events
  stack.X=X
  for(j in 2:d){
    stack.X=rbind(stack.X,X)
  }
  stack.W=W
  rows=rep(1:nrow(stack.W),cbind(1:nrow(stack.W),n)[,2])
  stack.W=stack.W[rows,]
  
  ### Initialize parameters mixture weights
  
  alpha.chain=matrix(0,maxT,K) ## canvas
  alpha.prior=rep(1,K)
  alpha=rdirichlet(1,alpha.prior)
  
  ## Initialize allocation variables z
  z=matrix(0,n,K) ## original allocation vectors
  index=rep(0,n)
  for (i in 1:n){
    index[i]=sample(1:K,1,prob=alpha)
    z[i,index[i]]=1
  }
  
  ## Build design matrix
  ## Compute linear predictor and probs
  ## initialize regr.coeff from prior
  design=array(0,c(K,n*d,L+Q))
  lin.pred=matrix(0,n*d,K)
  prior.m.coeff=matrix(0,L+Q,K)
  prior.var.coeff=array(1,c(L+Q,L+Q,K))
  beta.chain=array(0,c(K,L,maxT))
  gamma.chain=array(0,c(K,Q,maxT))
  temp.coeff=matrix(0,K,L+Q)
  for(k in 1:K){
    prior.m.coeff[,k]=rep(0,L+Q) 
    prior.var.coeff[,,k]=prior.var*diag(L+Q)
    temp.coeff[k,] <- rmvnorm(1,mean=prior.m.coeff[,k],sigma=prior.var.coeff[,,k])
    design[k,,]=cbind(stack.X,stack.W)
    lin.pred[,k]=design[k,,]%*%temp.coeff[k,]
  }
  probs=pnorm(lin.pred)
  probs[is.na(probs)]=0.00001
  probs[probs==0]=0.00001
  probs[probs==1]=0.99999
  
  ### canvas for post.prob of allocation
  prob.z.chain=array(0,c(maxT,n,K))
  z.chain=array(0,c(maxT,n,K))
  
  ### canvas for log.likelihood
  flag.na.z=log.lik=rep(0,maxT)
  
  ##
  out.coeff=matrix(0,K,L+Q)
  
  ##### Start MCMC
  for (tt in 1:maxT){
    
    ## Step 1. Sample z
    prob.z=matrix(0,n,K)
    for(k in 1:K){
      f.y.z=stack.y*log(probs[,k])+(1-stack.y)*log(1-probs[,k])
      prob.z[,k]=log(alpha[k])+rowSums(matrix(f.y.z,n,d,byrow=FALSE))
    }
    prob.z=exp(prob.z)
    prob.z=prob.z/rowSums(prob.z)
    if(any(is.na(prob.z))==1) flag.na.z[tt]=1
    prob.z=ifelse(is.na(prob.z),1/K,prob.z)
    prob.z.chain[tt,,]=prob.z
    z=t(apply(prob.z,1,function(x) {rmultinom(1,1,prob=x)}))
    index=apply(z,1,which.max)
    n_sums=colSums(z)
    # index=true.params$index
    
    full.index=cbind(seq_along(rep(index,d)),rep(index,d))
    
    ## Step 2. Retrieve n_star and Sample alpha_star
    alpha=rdirichlet(1,n_sums+alpha.prior)
    alpha.chain[tt,]=alpha
    
    
    ## Step 3. Sample regr.coeff; compute lin.pred and probs
    for(k in 1:K){
      if(n_sums[k]>0){
        out.coeff=probit_sampler(y=stack.y[full.index[,2]==k],
                             cov.mat=design[k,full.index[,2]==k,],
                             buff.coeff=temp.coeff[k,],
                             prior.m.coeff=prior.m.coeff[,k],
                             prior.var.coeff=prior.var.coeff[,,k])
        temp.coeff[k,]=out.coeff
      }
      beta.chain[k,,tt]=temp.coeff[k,1:L]
      gamma.chain[k,,tt]=temp.coeff[k,(L+1):length(temp.coeff[k,])]
    }
    
    ## Step 5. update preds and pnorms
    for(k in 1:K){
      lin.pred[,k]=design[k,,]%*%temp.coeff[k,]
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
  
  out=(list(prob.z.chain=prob.z.chain,
            alpha.chain=alpha.chain,
            beta.chain=beta.chain,
            gamma.chain=gamma.chain,
            last.cl=index,
            log.lik=log.lik,
            flag.na.z=flag.na.z))
}

########### Functions
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
