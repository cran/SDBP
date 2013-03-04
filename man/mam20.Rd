\name{mam20}
\docType{data}
\alias{mam20}
\title{Site-wise loglikelihoods scores for each of the mammalian mitochondrial 
amino acid sequences and the DNA sequences 12S and 16S rRNA genes}

\description{
Site-wise loglikelihoods scores for each of the 20 mammalian species' mitochondrial 
amino acid sequences (12 individual proteins) and the DNA sequences 12S and 
16S rRNA genes for 15 trees.}
\usage{data(mam20)}
\format{matrix of size \eqn{5879 \times 15}.}
\details{
  See \code{source} section for \code{data files}.
}

\examples{
\dontrun{
## Reading the data
library(SDBP)
data(mam20)### data named mam20 stored in the file mam20.rda

## Speedy double bootstrap probabilities for 15 trees
mam20.pv <- sdbp.default(mam20)
mam20.pv
}
}
\source{\url{http://www.bi.cs.titech.ac.jp/sdbp/}}
\references{
  Ren, A. et al. (2013) 
   "Assessing statistical reliability of phylogenetic trees via a speedy double bootstrap method",
  \emph{Molecular Phylogenetics and Evolution}, 
   67, 429-435.
}
\keyword{datasets}


