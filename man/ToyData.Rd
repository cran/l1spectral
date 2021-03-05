\name{ToyData}
\alias{ToyData}
\docType{data}
\title{
  Toy data for running the l1-spectral clustering algorithm
}
\description{
  An example of data for running the l1-spectral clustering algorithm.
}
\usage{ToyData}
\format{
  A list of three variables containing the adajcency matrix \code{A} of a 5-nodes graph, the adjacency matrix \code{A_hat} of a perturbed version of the same graph and the length of the two inherent clusters.
}
\value{
  No value returned, as this is a dataset.
}
\examples{
data(ToyData)
A <- ToyData$A
A_hat <- ToyData$A_hat
clusters <- ToyData$clusters
}
