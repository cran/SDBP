\name{sdbpk}
\alias{sdbpk}
\alias{print.sdbpk}
\title{Speedy double bootstrap probability for a specific phylogenetic tree}
\description{Calculate the speedy double bootstrap probability for a specific phylogenetic tree}
\usage{
sdbpk(dat,k,nb=10000,...)

\method{print}{sdbpk}(x,...)
}
\arguments{
  \item{dat}{a site-wise maximum log-likelihood matrix.
             Each column is for a tree, each row is for a site-wise maximum log-likelihood.}
  \item{k}{a specific tree number in the .mt file.}
  \item{nb}{the number of first-tier replicates is sampled from the multivariable normal
           distribution
           \eqn{N_T(\hat \mu,\Sigma)}. 
           T is the number of trees,
           \eqn{\hat \mu} is the T-degree vector,
           it is a projection of the maximum log-likelihood vector, which consists of
           the maximum log-likelihood for each tree,
           \eqn{\Sigma} is the unbiased covariance matrix of the maximum log-likelihood vector.}
  \item{x}{an object of class \code{sdbpk}.} 
  \item{\dots}{not used.}
}

\value{
  An object of class \code{sdbpk}, basically a list including elements
  \item{value}{a named numeric of the speedy double bootstrap probability for a tree k.}
  \code{print} method for an object of class \code{sdbpk}; prints the speedy double bootstrap probability
   for a tree k.}
\author{Ren Aizhen}
\seealso{\code{\link{sdbp.default}}.}
\examples{
library(SDBP)
data(mam20)### data named mam20 stored in the file mam20.rda
### 1 is the tree number ordered in the file trees15.tpl
result<-sdbpk(mam20,1)
result

### the output 
## Call:
## sdbpk(dat = data, k = 1)

##    t1 
## 0.7489   

}

\keyword{speedy double bootstrap probability for a specific phylogenetic tree}







