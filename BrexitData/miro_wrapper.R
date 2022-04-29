source('miro_int.R')
source('miro_NA.R')
#### Wrapper function for miro
miro<-function(y,X,W,K,prior.var=5,maxT=5000,seed=42,
                  nullClust=FALSE){
  if(is.data.frame(y)){
    y <- as.matrix(y)
  }
  if(is.null(X)){
    if(is.null(W)){
      # no covariates, call miro_A (actors) with intercept only
      X <- as.matrix(rep(1,dim(y)[1]))
      source('miro_actor.R')
      miro_A(y,X,K,prior.var,maxT,seed,nullClust)
    }else{
      # event covariates only, call miro_E
      source('miro_event.R')
      miro_E(y,W,K,prior.var,maxT,seed,nullClust)
    }
  }else{
    if(is.null(W)){
      # actor covariates only, call miro_A
      source('miro_actor.R')
      miro_A(y,X,K,prior.var,maxT,seed,nullClust)
    }else{
      # actor and event covariates, call miro_AE
      source('miro_actor_event.R')
      miro_AE(y,X,W,K,prior.var,maxT,seed,nullClust)
    }
  }
}