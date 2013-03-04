\name{bp}
\alias{bpk}
\alias{bp}
\alias{print.bp}
\title{Bootstrap probabilities for phylogenetic trees}
\description{Calculate the bootstrap probabilities for phylogenetic trees}
\usage{
bpk(dat,k,nb=10000,...)
bp(dat,nb=10000,...)

\method{print}{bp}(x,...)
}
\arguments{
  \item{dat}{a site-wise maximum log-likelihood matrix.
             Each column is for a tree, each row is for a site-wise maximum log-likelihood.}
  \item{k}{a specific tree number in the .mt file.}
  \item{nb}{the number of regular bootstrap replicates sampled from a multivariable normal
           distribution
           \eqn{N_T(y,\Sigma)}. 
           T is the number of trees,
           y is the T-degree vector,
           which is the maximum log-likelihood vector consisting of
           the maximum log-likelihood for each tree,
           \eqn{\Sigma} is the unbiased covariance matrix of the maximum log-likelihood vector.}
  \item{\dots}{not used.}
   \item{x}{an object of class \code{bp}.} 
}

\value{
  An object of class \code{bpk}, basically a list including elements 
  \item{BP}{the bootstrap probability for a specific tree k.}
  \item{t}{a specific tree number k.}
  An object of class \code{bp}, basically a list including element
  \item{pvs}{the bootstrap probabilities for the trees.}
  \code{print} method for an object of class \code{bp}; prints the bootstrap probabilities
   for the trees.}
  
\author{Ren Aizhen}
\examples{
library(SDBP)
data(mam20) ### data named mam20 stored in the file mam20.rda
result1<-bpk(mam20,2)
result1
\dontrun{
result2<-bp(mam20)
result2
}
### the output of result1
## $BP
## [1] 0.0292

## $t
## [1] 2

### the output of result2 is in decreasing order of log-likelihood for 
###   each tree  
## Call:
## bp(dat = data)

## Bootstrap probabilities:
##     t1     t4     t3     t7     t2     t5     t9    t11    t10      
## 0.4735 0.1934 0.1197 0.0927 0.0268 0.0422 0.0011 0.0256 0.0183 
##         t6     t8    t14    t13    t12    t15
##  0.0010 0.0003 0.0000 0.0000 0.0000 0.0000 
}

\keyword{bootstrap probabilities for the phylogenetic trees}







