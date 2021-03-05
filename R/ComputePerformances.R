#' @title Compute the performances of the l1-spectral clustering algorithm
#' @description This function computes the performances of the l1-spectral clustering algorithm in terms of Normalized Mutualized Information (NMI).
#' @param Results Output of the function \code{l1_spectralclustering()}.
#' @param A The adjacency matrix of the graph to cluster.
#' @importFrom NMI NMI
#' @seealso \code{\link{l1_spectralclustering}}, \code{\link{l1spectral}}.
#' @return The Normalized Mutualized Information (NMI) score.
#' @author Camille Champion, Magali Champion
#' @export
#' @examples
#'  #############################################################
#'  # Generating toy data
#'  #############################################################
#'  Data <- CreateDataSet(k=3, n=20, p=list(p_inside=0.1,p_outside=0.1))
#'

ComputePerformances <- function(Results, A){
  # Results: output of the function l1_spectralclustering()
  # A: true adjacency matrix

  # Output: NMI score

  # first, find the clusters in the adjacency matrix
  graph <- graph_from_adjacency_matrix(A,mode="undirected")
  clusters <- components(graph)$membership
  clusters <- cbind(c(1:length(clusters)),clusters)

  if (!is.null(ncol(Results$comm))){
    clus_est <- Results$comm%*%c(1:ncol(Results$comm))
  } else {
    clus_est <- Results$comm
  }
  clus_est <- cbind(c(1:length(clus_est)),clus_est)

  score <- NMI(clus_est,clusters)
  score <- score$value
  return(score)
}
