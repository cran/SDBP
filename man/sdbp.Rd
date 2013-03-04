\name{sdbp}
\alias{sdbp.default}
\alias{print.sdbp}
\alias{summary.sdbp}
\title{Speedy double bootstrap probabilities for candidate phylogenetic trees}
\description{Calculates the speedy double bootstrap probabilities for candidate trees}
\usage{

\method{sdbp}{default}(dat,nb=10000,...)

\method{print}{sdbp}(x,...)	

\method{summary}{sdbp}(object,...)

}
\arguments{
  \item{dat}{a site-wise maximum log-likelihood matrix. Each column is for a tree. 
    Each row is for a site-wise maximum log-likelihood.}
  \item{nb}{the number of first-tier replicates sampled from a multivariable normal
           distribution
           \eqn{N_T(\hat \mu,\Sigma)}. 
           T is the number of trees,
           \eqn{\hat \mu} is the T-degree vector,
           it is a projection of the maximum log-likelihood vector, which consists of
           the maximum log-likelihood for each tree,
           \eqn{\Sigma} is the unbiased covariance matrix of the maximum log-likelihood vector.}
   \item{x}{an object of class \code{sdbp}.}
  \item{object}{an object of class \code{sdbp}.}  
  \item{\dots}{not used.}
}
\value{
  class \code{sdbp} returns an object, basically a list including elements
  \item{pvs}{a named numeric of speedy double bootstrap probabilities for all candidate trees.}
  \item{nb}{same as argument \code{nb}.}
  \code{print} method for an object of class \code{sdbp}; prints the speedy double bootstrap probabilities
   for the trees.
  \code{summary} method for an object of class \code{sdbp}; gives the speedy double bootstrap probabilities
   and their standard errors for the trees.
}
\author{Ren Aizhen}
\seealso{\code{\link{sdbpk}}.}
\examples{
\dontrun{
library(SDBP)
data(mam20)### data named mam20 stored in the file mam20.rda
result <- sdbp.default(mam20)
result
summary(result)
}

### The results are output in decreasing order of log-likelihood for each tree
## > result
## Call:
## sdbp.default(dat = data)

## Speedy double bootstrap probabilities:
##    t1     t4     t3     t7     t2     t5     t9    t11    t10     t6     t8 
## 0.7575 0.4258 0.3762 0.3330 0.2984 0.2592 0.1638 0.2164 0.1801 0.0472 0.0346 
##    t14    t13    t12    t15 
## 0.0101 0.0107 0.0121 0.0009 

## > summary(result)
## Call:
## sdbp.default(dat = data)

##     StdErr   sDBP    
## t1  0.0043 0.7575    
## t4  0.0049 0.4258    
## t3  0.0048 0.3762    
## t7  0.0047 0.3330    
## t2  0.0046 0.2984    
## t5  0.0044 0.2592    
## t9  0.0037 0.1638    
## t11 0.0041 0.2164    
## t10 0.0038 0.1801    
## t6  0.0021 0.0472 *  
## t8  0.0018 0.0346 *  
## t14 0.0010 0.0101 *  
## t13 0.0010 0.0107 *  
## t12 0.0011 0.0121 *  
## t15 0.0003 0.0009 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
   
}
\keyword{speedy double bootstrap probabilities for the phylogenetic trees}







