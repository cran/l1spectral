\name{l1spectral-package}
\alias{l1spectral}
\docType{package}
\title{
  Description of the package
}
\seealso{\code{\link{l1_spectralclustering}}}

\description{
  \packageDescription{l1spectral}
}
\details{
l1-spectral clustering is an l1-penalized version of the spectral clustering algorithm, which aims at robustly detecting cluster structure of perturbed graphs by promoting sparse eigenbases solutions of specific l1-minimization problems.

  The DESCRIPTION file:
    \packageDESCRIPTION{l1spectral}
}
\author{
  \packageAuthor{l1spectral}
}
\keyword{Unsupervised learning - Spectral clusterin - l1-penalty - Biological networks}
\references{
C. Champion, M. Champion, M. Blazère, R. Burcelin, J.M. Loubes, l1-spectral clustering algorithm: a robust spectral clustering using Lasso regularization, Preprint (2021).
}
\examples{
 #####################################################
 # Performing the l1-spectral clustering on the graph
 #####################################################

 data(ToyData)

 # if desired, the number of clusters and representative elements can be provided,
 # otherwise remove
 results2 <- l1_spectralclustering(A = ToyData$A_hat, pen = "lasso")
 results2$comm

 # when desired, the number of clusters and representative elements can also be provided
 \donttest{results2 <- l1_spectralclustering(A = ToyData$A_hat, pen = "lasso",
              k=2, elements = c(1,4))}
}
