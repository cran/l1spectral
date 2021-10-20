#' @importFrom stats lm
#' @importFrom grDevices rainbow
#' @importFrom Rcpp evalCpp
#' @useDynLib l1spectral, .registration = TRUE
#'
#' @title Run the l1-spectral clustering algorithm
#' @description This function runs the l1-spectral algorithm, an l1-penalized version of the spectral clustering that aims at robustly clustering perturbed graphs.
#'
#' @param A The adjacency matrix of the graph to cluster.
#' @param k True number of clusters (not necessarily needed). If not provided, k is chosen by spectral eigengap.
#' @param elements The representative elements of the clusters (not necessary needed). If not provided, index are chosen using the betweeness centrality score.
#' @param pen The penalty (to be chosen among "lasso" or "thresholdedLS").
#' @param stab TRUE/FALSE indicated whether the indices should be stabilized (TRUE by default)
#'
#' @seealso \code{\link{ComputePerformances}}, \code{\link{l1spectral}}.
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{comm}}{ The community matrix,}
#' \item{\code{structure}}{ The structure of the graph to cluster,}
#' \item{\code{clusters}}{ The number of clusters,}
#' \item{\code{elements}}{ The chosen representative elements of the clusters.}
#' }
#' @export
#'
#' @author Camille Champion, Magali Champion
#' @examples
#'  #####################################################
#'  # Performing the l1-spectral clustering on the graph
#'  #####################################################
#'
#'  data(ToyData)
#'
#'  # if desired, the number of clusters and representative elements can be provided, otherwise, remove
#'  results2 <- l1_spectralclustering(A = ToyData$A_hat, pen = "lasso")
#'  results2$comm
#'
#'  # when desired, the number of clusters and representative elements can also be provided
#'  \donttest{results2 <- l1_spectralclustering(A = ToyData$A_hat, pen = "lasso",
#'              k=2, elements = c(1,4))}

l1_spectralclustering <- function(A, k = NULL, elements = NULL, pen, stab = TRUE){
  # A: the matrix we aim at clustering (adjacency matrix, e.g. from CreateDataSet())
  # k: true number of clusters (not necessary needed). If not precised, k is chosen by spectral eigengap
  # elements: representative elements of the clusters (not necessary needed). If not precised, index are chosen using the betweeness centrality score.
  # pen: penalty (to be chosen among "lasso" or "thresholdedLS")
  # stab: TRUE/FALSE, should the indices be stabilized?

  # outputs:
  #   - comm: the community matrix
  #   - structure: structure of the graph
  #   - clusters: number of clusters
  #   - elements: the representative elements of the clusters

  # 1st step: finding the connected components
  Structure <- FindStructure(A)

  # 2nd step: finding the optimal number of clusters (only if k is not provided)
  clusters <- FindNbrClusters(A, structure  = Structure, k = k)

  # 3rd step: finding the representative elements of the clusters
  Elements <- FindElement(A = A, structure = Structure, clusters = clusters, elements = elements)

  # 4th step: running the l1-spectral clustering algorithm on each connected component (each of them are treated independtly)
  comm <- matrix(0,nrow=ncol(A),ncol=clusters$nbr_clusters_total)
  S <- cumsum(unlist(clusters$nbr_clusters)[-length(clusters$nbr_clusters)])
  for (i in (1:length(Structure$groups))){
    Atmp <- A[Structure$groups[[i]],Structure$groups[[i]]]
    clusters_tmp <- clusters$nbr_clusters[[i]]
    indices_tmp <- Elements$indices[[i]][which(Elements$indices[[i]]%in%Structure$groups[[i]])]
    indices_tmp <- match(indices_tmp,Structure$groups[[i]])
    score_tmp <- Elements$score[[i]]
    names(score_tmp) <- paste0("Node",match(as.numeric(substring(names(Elements$score[[i]]),5)),Structure$groups[[i]]))
    Elements_tmp <- list(score = score_tmp,indices = indices_tmp)

    results <- l1_spectral(A = Atmp, k = clusters_tmp, elements = Elements_tmp, pen = pen, stab=stab)

    print(paste0("Component ",i," clustered."))

    if (!is.null(ncol(results))){
      if (ncol(results)!=clusters$nbr_clusters[[i]]){
        clusters$nbr_clusters[[i]] <- ncol(results)
        clusters$nbr_clusters_total <- sum(unlist(clusters$nbr_clusters))
      }
      S <- cumsum(unlist(clusters$nbr_clusters))
    } else {
      if (clusters$nbr_clusters[[i]]!=1){
        clusters$nbr_clusters[[i]] <- 1
        clusters$nbr_clusters_total <- sum(unlist(clusters$nbr_clusters))
      }
      S <- cumsum(unlist(clusters$nbr_clusters))
    }
    if (i==1){
      comm[Structure$groups[[i]],1:S[1]] <- results
    } else {
      comm[Structure$groups[[i]],(S[i-1]+1):S[i]] <- results
    }
  }

  if (length(which(colSums(comm)==0))>0){
    comm <- comm[,-which(colSums(comm)==0)]
  }
  comm[comm>0] <- 1
  return(list(comm=comm,Structure=Structure,clusters=clusters,Elements =Elements))
}
