#' @importFrom igraph betweenness
#'
#' @title Find the representative elements of the clusters
#'
#' @description This internal function of the l1-spectral clustering algorithm finds representative elements of the clusters, that is nodes belonging to the clusters.
#' @param A The adjacency matrix
#' @param structure Output of the function \code{FindStructure()}.
#' @param clusters Output of the function \code{FindNbrClusters()}.
#' @param elements The representative elements of the clusters (not necessary needed). If not provided, chosen using the betweeness centrality score.
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{score}}{ The edge betweenness score of all nodes,}
#' \item{\code{Nodes}}{ Vector of the representative elements.}
#' }
#' @export
#'
#' @seealso \code{\link{l1_spectralclustering}}, \code{\link{l1spectral}}.
#' @author Camille Champion, Magali Champion
#' @examples
#'  ######################################################
#'  # Finding the representative elements of the clusters
#'  ######################################################
#'
#'  # 1st: create data (not perturbed graph)
#'  Data <- CreateDataSet(k=3, n=20, p=list(p_inside=0,p_outside=0))
#'
#'  # 2nd: find the structure of the graph
#'  Structure <- FindStructure(Data$A_hat)
#'
#'  # 3rd: find the optimal number of clusters (here, 3 clusters)
#'  Clusters <- FindNbrClusters(A = Data$A_hat, structure = Structure, k=3)
#'
#'  # 4th: find the representative elements of the clusters
#'  Elements <- FindElement(A = Data$A_hat, structure = Structure, clusters = Clusters)
#'  # if elements is not provided, the representative elements of each component are chosen
#'  # by maximizing the edge betweenness score
#'
#'  Elements <- FindElement(A = Data$A_hat, structure = Structure,
#'                    clusters = Clusters, elements = c(1,5,12))

FindElement <- function(A, structure, clusters, elements = NULL){
  # A: the adjacency matrix of the graph
  # structure: Structure: already existing connected components
  # clusters: output of the function FindNbrClusters()
  # elements: representative elements of the clusters (not needed)
  # Outputs:
  #   - score: the edge betweenness score of all nodes
  #   - Nodes: vector of the representative elements

  n <- ncol(A)

  between <- function(comm){
    # betweeness centrality score on subgraphs of a graph
    graph_small <- graph_from_adjacency_matrix(A[comm,comm],mode="undirected")
    b <- betweenness(graph_small,directed = FALSE,normalized=TRUE)
    names(b) <- paste0("Node",comm)
    I <- order(b,decreasing = TRUE)
    b <- b[I]
    return(b)
  }

  betweenness <- lapply(structure$groups,between)

  if (is.null(elements)){
    # case 1: no representative elements

    Nodes <- c()
    for (i in (1:length(betweenness))){
      b <- betweenness[[i]][1:unlist(clusters$nbr_clusters)[i]]
      b <- as.numeric(substring(names(rev(b)), 5))
      Nodes <- c(Nodes,list(b))
    }
    names(Nodes) <- names(structure$groups)

  } else {
    # case 2: provided representative elements
    if (length(elements)<clusters$nbr_clusters_total){
      print("The number of representative elements does not coincide with the number of clusters. Please add some elements to re-run the code. By default, new elements will be chosen.")

      Nodes <- c()
      for (i in (1:length(betweenness))){
        b <- betweenness[[i]][1:unlist(clusters$nbr_clusters)[i]]
        b <- as.numeric(substring(names(rev(b)), 5))
        Nodes <- c(Nodes,list(b))
      }
      names(Nodes) <- names(structure$groups)
    } else{
      test <- lapply(structure$groups,function(x){
        length(which(elements%in%x))
      })
      test2 <- unlist(test) - unlist(clusters$nbr_clusters)

      if (length(which(test2<0))>0){
        print("At least one of the connected components does not have enough representative element. New elements will be chosen.")
        Nodes <- c()
        for (i in (1:length(betweenness))){
          b <- betweenness[[i]][1:unlist(clusters$nbr_clusters)[i]]
          b <- as.numeric(substring(names(rev(b)), 5))
          Nodes <- c(Nodes,list(b))
        }
        names(Nodes) <- names(structure$groups)
      } else {
        Nodes <- c()
        for (i in (1:length(betweenness))){
          index_gr <- elements[which(elements %in% structure$groups[[i]])]
          b <- betweenness[[i]][paste0("Node",index_gr)]
          if (length(b)>1){
            b <- sort(b,decreasing = TRUE)
          }
          if (test2[i] > 0){
            b <- b[1:clusters$nbr_clusters[[i]]]
          }
          b <- as.numeric(substring(names(rev(b)), 5))
          Nodes <- c(Nodes,list(b))
        }
        names(Nodes) <- names(structure$groups)
      }
    }
  }
  return(list(score=betweenness,indices=Nodes))
}
