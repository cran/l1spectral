#' @importFrom igraph components
#' @importFrom igraph groups
#' @title Find the structure of the graph from the adjacency matrix
#'
#' @description This internal function of the spectral clustering algorithm finds the structure of the graph to cluster (number of nodes and connected components).
#' @param A The adjacency matrix
#'
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{graph}}{ igraph object derived from A,}
#' \item{\code{groups}}{ List of connected components and corresponding nodes.}
#' }
#' @seealso \code{\link{l1_spectralclustering}}, \code{\link{l1spectral}}.
#' @export
#' @author Camille Champion, Magali Champion
#' @examples
#'  ###############################################################
#'  # Finding the structure of the graph from the adjacency matrix
#'  ###############################################################
#'
#'  # 1st example: non-perturbed graph
#'  Data <- CreateDataSet(k=3, n=20, p=list(p_inside=0,p_outside=0))
#'
#'  Structure <- FindStructure(Data$A_hat)
#'  Structure$groups # the graph is not perturbed, there are 3 connected components
#'
#'  # 2nd example: highly-perturbed graph
#'  Data <- CreateDataSet(k=3, n=20, p=list(p_inside=0.5,p_outside=0.5))
#'
#'  Structure <- FindStructure(Data$A_hat)
#'  Structure$groups # the graph is higlhy perturbed, there are less than 3 connected components


FindStructure <- function(A){
  # A: adjacency matrix of the graph
  # Output: list of two elements
  #   - graph: igraph object with the graph
  #   - groups: list of connected components and corresponding nodes

  # define the graph
  graph <- graph_from_adjacency_matrix(A,mode="undirected")

  # find the connected components
  clu <- components(graph)
  groups <- groups(clu)
  if (length(groups)>1){
    groups <- groups[order(sapply(groups,length))]
  } else {
    groups <- groups
  }
  nbr_comp <- length(groups)
  names(groups) <- paste0("Component",c(1:nbr_comp))

  if (nbr_comp==1){
    print("The graph has only one connected component.")
  } else {
    print(paste0("The graph have ",nbr_comp," connected components. Each of them will be clustered using the l1-spectral clustering."))
  }
  return(list(graph=graph,groups=groups))
}
