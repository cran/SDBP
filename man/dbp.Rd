\name{dbpk}
\alias{dbpk}
\title{Double bootstrap probability for a specific phylogenetic tree}
\description{Calculate the double bootstrap probability for a specific phylogenetic tree}
\usage{

dbpk(dat,k,nb1=1000,nb2=1000,nb=5000,...)
}
\arguments{
  \item{dat}{a site-wise maximum log-likelihood matrix.
             Each column is for a tree, each row is for a site-wise maximum log-likelihood.}
  \item{k}{a specific tree number in the .mt file.}
  \item{nb1}{the number of first-tier bootstrap replicates sampled 
             from a multivariable normal distribution
           \eqn{N_T(\hat\mu,\Sigma)}. 
           T is the number of trees,
           \eqn{\hat\mu} is the T-degree vector,
           which is a projection of the maximum log-likelihood vector, consisting of the
           maximum log-likelihood for each tree,
           \eqn{\Sigma} is the unbiased variance matrix of the maximum log-likelihood vector.}
  \item{nb2}{the number of second-tier bootstrap replicates sampled
             from a multivariable normal distribution
             \eqn{N_T(y^*,\Sigma)}. 
             T is the number of trees,
             \eqn{y^*} is the T-degree vector,
             which are first-tier bootstrap replicates,
            \eqn{\Sigma} is the unbiased variance matrix of the maximum log-likelihood vector.}
  \item{nb}{the number of regular bootstrap replicates sampled from a multivariable 
            normal distribution
           \eqn{N_T(y,\Sigma)}. 
           T is the number of trees,
           y is the T-degree vector,
           which is the maximum log-likelihood vector, consisting of the 
           maximum log-likelihood for each tree,
           \eqn{\Sigma} is the unbiased variance matrix of the maximum log-likelihood vector.}
\item{\dots}{not used.}
}

\value{
  An object of class \code{dbpk}, basically a list including elements
  \item{DBP}{the double bootstrap probability for a tree k.}
  \item{t}{a specific phylogenetic tree numbered k.}
 }
\author{Ren Aizhen}
\examples{
\dontrun{
library(SDBP)
data(mam20)### data named mam20 stored in the file mam20.rda
### 3 is the tree number ordered in file trees15.tpl 
result<-dbpk(mam20,3)
result
}
### the output of result 

## $DBP
##    t3 
## 0.317 

## $t
## [1] 3
}

\keyword{double bootstrap probability for a specific phylogenetic tree}







