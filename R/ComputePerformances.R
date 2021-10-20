#' @title Compute the performances of the l1-spectral clustering algorithm
#' @description This function computes the performances of the l1-spectral clustering algorithm in terms of Normalized Mutualized Information (NMI).
#' @param Results Output of the function \code{l1_spectralclustering()}.
#' @param A The adjacency matrix of the graph to cluster.
#' @importFrom aricode NMI
#' @importFrom aricode AMI
#' @importFrom aricode ARI
#' @seealso \code{\link{l1_spectralclustering}}, \code{\link{l1spectral}}.
#' @return The Normalized Mutualized Information (NMI), Adjusted Mutualized Information (AMI) and Adjusted Rand Index (ARI) scores.
#' @author Camille Champion, Magali Champion
#' @export
#' @examples
#'  #############################################################
#'  # Computing the performances
#'  #############################################################
#'
#'  data(ToyData)
#'
#'  results <- l1_spectralclustering(A = ToyData$A_hat, pen = "lasso",
#'              k=2, elements = c(1,4))
#'
#'  ComputePerformances(Results=results,A=ToyData$A)
#'

ComputePerformances <- function(Results, A){
  # Results: output of the function l1_spectralclustering()
  # A: true adjacency matrix

  # Output: NMI, AMI and ARI scores

  # first, find the clusters in the adjacency matrix
  graph <- graph_from_adjacency_matrix(A,mode="undirected")
  clusters <- components(graph)$membership

  if (!is.null(ncol(Results$comm))){
    clus_est <- Results$comm%*%c(1:ncol(Results$comm))
  } else {
    clus_est <- Results$comm
  }
  clus_est <- as.vector(clus_est)

  NMI <- NMI(clus_est,clusters)
  AMI <- AMI(clus_est,clusters)
  ARI <- ARI(clus_est,clusters)

  return(score=list(NMI=NMI,AMI=AMI,ARI=ARI))
}
