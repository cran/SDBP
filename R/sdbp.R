##  SDBP: R package for speedy double bootstrap
##  Copyright (C) 2013 Aizhen Ren
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation.
##
######################################################################
### MAIN: TESTING FOR PHYLOGENETIC TREE SELECTION PROBLEM via SPEEDY DOUBLE BOOTSTRAP METHOD
## dat: data matrix (produced by read.mt)
## nb: numbers of replicates
## mypava: calculate the projection for  maximum log-likelihood vector y using pava method
## k: tree number in the .mt file  
### 
### calculate the projection of y for hypo H_k which is tree k is best tree 
##
### create h dimensions normal random variables
## mu : mean vector 
## A : unbiased covariance matrix 
## r : number of normal randam variable
## h : dimensions of the normal random variables 
rknorm <- function(mu, A, r, h) {
    	U <- svd(A)$u
    	V <- svd(A)$v
    	D <- diag(sqrt(svd(A)$d))
    	B <- U %*% D %*% t(V) # sqrt of matrix A
    	w <- c()
    	for (i in 1:r)
    	w <- append(w, list(mu + B%*%cbind(rnorm(h))))
    	return(w)
    }
### calculate the projection of vector y
## y : log-likelihood vector 
## k : tree's number  
mypava <- function(y,k){
    count <- 1 
    for(i in seq(along=y)){
       if(y[k]<y[i]){
          y[k] <- (count*y[k]+y[i])/(count+1)
          count <- count+1
       }
    }
for(i in seq(along=y)[-k]){
   y[i] <- min(c(y[k],y[i]))
   }
   return(y)
   }
### for singel tree
sdbpk <- function(dat,k,nb = 10000,...) 
  { y <- apply(dat,2,sum) # log-likelihood vector
    nsize <- nrow(dat) # sample size of dataset
    covar <- nsize*var(dat) # unbiased covariance matrix of vector y  
    h <- ncol(dat) # number of trees  
    names(y) <- NULL
    d <- numeric(h)  
    for(i in 1:h){
      d[i] <- max(y[-i])-y[i]#the signed distance from original log-likelihood vector y for tree i
    }
    mybyy_2 <- function(haty,n,matr){### mybyy_2 are replications of haty for tree k
 	  by <- matrix(0,n,h)
        for(i in 1:n) {
          by[i,] <- rknorm(haty,matr,1,h)[[1]]#replicates of haty for tree k
        }
        return(by)
    }
    hatmu <- mypava(y,k)#y's projection for tree k by using pava algorithm
    by <- mybyy_2(hatmu,nb,covar)#first-tier bootstrap replicates
    d_as <- numeric(nb)## the signed distance of replications by
    for(i in 1:nb){
       d_as[i] <- max(by[i,][-k])-by[i,][k]
    }
    sDBP <- sum(d_as > d[k])/nb
    sdbpk <- list("value"=sDBP)
    names(sdbpk$value) <- paste("t",k,sep="")
    sdbpk$call <- match.call()
    class(sdbpk) <- "sdbpk" 
    sdbpk
}
print.sdbpk <- function(x,...){
     	cat("Call:\n")
    	print(x$call)
    	cat("\n")
    	print(x$value)
}

sdbp.default <- function(dat,nb = 10000,...){
    	 if(ncol(dat) < 2) stop("should be more than two trees")
    	 y <- apply(dat,2,sum)###log-likelihood vector
	 originalt <- rev(order(y))
       dat <- dat[,originalt]###reorder the trees in decreasing order of log-likelihood for each tree
       h <- ncol(dat)
       y <- apply(dat,2,sum)###original maximum log-likelihood of each tree
       fun <- function(dat)sapply(1:h,function(x)sdbpk(dat=dat,x,nb=nb,...))
       est <- fun(dat)
       est$pvs <- as.numeric(est[1,])
       names(est$pvs) <- paste("t",originalt,sep="")
       est$nb <- nb
       est$call <- match.call()
       class(est) <- "sdbp" 
       est
    }
print.sdbp <- function(x,...){
      cat("Call:\n")
      print(x$call)
      cat("\nSpeedy double bootstrap probabilities:\n")
      print(x$pvs)
}
### summary.sdbp
summary.sdbp <- function(object,...){ 
     pval <- object$pvs
     nb <- object$nb
     se <- sqrt(pval*(1-pval)/nb)
     TAB <- cbind( StdErr= round(se,4),
                   sDBP=pval)
     res <- list(call=object$call,
                 coefficients=TAB)
     class(res)<-"summary.sdbp"
     res
}
print.summary.sdbp <- function(x,...){
      cat("Call:\n")
      print(x$call)
      cat("\n")
      printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)
      }
### calculate the bootstrap probability for tree k 
## dat: observerd data
## k: tree number in the .mt file
## nb: bootstrap replicates

 bpk <- function(dat,k,nb=10000,...){
     y <- apply(dat,2,sum)
     nsize <- nrow(dat)  # sample size of dataset
     h <- ncol(dat)
     covar <- nsize*var(dat) # Unbiased covariance matrix
     mu <- matrix(0,nb,h) # calculate bootstrap replicates
     for(i in 1:nb)
     {
       mu[i,] <- rknorm(y,covar,1,h)[[1]]
     }
     BP <- sum(apply(mu,1,which.max)==k)/nb
     list("BP"=BP,"t"=k)
}
### calculate the bootstrap probabilities for all tree
bp <- function(dat,nb = 10000,...){
     nsize <- nrow(dat)
     h <- ncol(dat)
     y <- apply(dat,2,sum)###original maximum log-likelihood of each tree
     originalt <- rev(order(y))
     dat <- dat[,originalt]
     fun <- function(dat)sapply(1:h,function(x)bpk(dat=dat,x,nb=nb,...))
     bp <- fun(dat)
     bp$pvs <- as.numeric(bp[1,])
     names(bp$pvs) <- paste("t",originalt,sep="")
     bp$call <- match.call()
     class(bp) <- "bp" 
     bp
}
print.bp <- function(x,...){
    cat("Call:\n")
    print(x$call)
    cat("\nBootstrap probabilities:\n")
    print(x$pvs)
   }
### calculate the double bootstrap probability for tree k 
## dat: observerd data
## k: tree number in the .mt file
## nb1: first-tier bootstrap replicates
## nb2: seconde-tier bootstrap replicates
## nb: seconde-tier bootstrap replicates for BP
dbpk <- function(dat,k,nb1=1000,nb2=1000,nb=5000,...){  
    if(ncol(dat) < 2) stop("should be more than two trees")
    y <- apply(dat,2,sum) # redefinite the y 
    nsize <- nrow(dat) # sample size of dataset
    covar <- nsize*var(dat) # Unbiased Covariance matrix of vector y  
    h <- ncol(dat) # number of trees 
    names(y) <- NULL
    mybyy_2 <- function(haty,n,covar){###mybyy_2 are replications of haty
        by <- matrix(0,n,h)
        for(i in 1:n) {
          by[i,] <- rknorm(haty, covar, 1, h)[[1]]#replicates of haty for tree k
        }
        return(by)
    }
    hatmu <- mypava(y,k)#y's projection for tree k by using pava algorithem
    by <- mybyy_2(hatmu,nb1,covar)#first-tier bootstrap replicates
    byy_y <- mybyy_2(y,nb,covar)#second-tier bootstrap replicates for y
    BP <- sum(apply(byy_y,1,which.max) == k)/nb
    componentDBP <- numeric(nb1)
    for(i in 1:nb1){
      byy <- mybyy_2(by[i,],nb2,covar)#second-tier bootstrap replicates for by[i]
      componentDBP[i] <- sum(apply(byy,1,which.max) == k)/nb2
    }
    DBP <- sum(componentDBP<= BP )/nb1
    names(DBP) <- paste("t",k,sep="")
    list("DBP"=DBP,"t"=k)
    }

  
   
     

   



   
  
   
   

    
   
  
 

