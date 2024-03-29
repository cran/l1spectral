% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ComputePerformances.R
\name{ComputePerformances}
\alias{ComputePerformances}
\title{Compute the performances of the l1-spectral clustering algorithm}
\usage{
ComputePerformances(Results, A)
}
\arguments{
\item{Results}{Output of the function \code{l1_spectralclustering()}.}

\item{A}{The adjacency matrix of the graph to cluster.}
}
\value{
The Normalized Mutualized Information (NMI), Adjusted Mutualized Information (AMI) and Adjusted Rand Index (ARI) scores.
}
\description{
This function computes the performances of the l1-spectral clustering algorithm in terms of Normalized Mutualized Information (NMI).
}
\examples{
 #############################################################
 # Computing the performances
 #############################################################

 data(ToyData)

 results <- l1_spectralclustering(A = ToyData$A_hat, pen = "lasso",
             k=2, elements = c(1,4))

 ComputePerformances(Results=results,A=ToyData$A)

}
\seealso{
\code{\link{l1_spectralclustering}}, \code{\link{l1spectral}}.
}
\author{
Camille Champion, Magali Champion
}
