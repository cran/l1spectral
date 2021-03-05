#' @importFrom graphics plot
#' @importFrom  graphics par
#' @importFrom stats rbinom
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom Matrix bdiag
#' @title Create data set
#' @description This function generates toy data that can be used to run the l1-spectal clustering algorithm: the adjacency matrix of a graph with \code{n} nodes and its perturbed version.
#' @param k True number of clusters.
#' @param n Number of nodes.
#' @param p List of probabilities of perturbations (inside and outside clusters).
#' @param print.plot TRUE/FALSE indicated whether the graph should be plotted.
#' @param ClustersLength Length of the \code{k} clusters (not necessary needed). If not provided, randomly chosen in such a way that \code{sum(ClustersLength)=n}.
#' @seealso \code{\link{l1_spectralclustering}}, \code{\link{l1spectral}}.
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{A}}{ Adjacency matrix of the generated graph.}
#' \item{\code{A_hat}}{ Adjacency matrix of the perturbed version of the generated graph.}
#' \item{\code{ClustersLength}}{ Length of the \code{k} clusters.}
#' }
#' @author Camille Champion, Magali Champion
#' @export
#' @examples
#'  #############################################################
#'  # Generating toy data
#'  #############################################################
#'  Data <- CreateDataSet(k=3, n=20, p=list(p_inside=0.1,p_outside=0.1))
#'
#'  # Data is a list of three objects:
#'  # - Data$A is an nxn matrix corresponding to the adjacency matrix of a graph
#'  # with n nodes and k clusters,
#'  # - Data$A_hat is a perturbed version of this graph with a probability
#'  # p_inside of removing an edge inside clusters and
#'  # p_outside of adding an edge between clusters,
#'  # - Data$ClustersLength is a vector indicating the length of the clusters.
#'
#'  Data <- CreateDataSet(k=3, n=20, p=list(p_inside=0.1,p_outside=0.1), print.plot=TRUE)
#'
#'  # The same as above but the true graph and its perturbed version are both plotted.

CreateDataSet <- function(k, n, p, print.plot=TRUE, ClustersLength = NULL){
  # k: true number of clusters
  # n: number of nodes
  # perturbation: list of perturbations, p1 (inside) and p2 (outside)
  # print.plot: TRUE/FALSE should the graph be plotted?
  # ClustersLength: length of the k clusters (not needed). If not precised, randomly chosen in such a way that sum(ClustersLength)=n

  # Outputs:
  # A: non-perturbed adjacency matrix
  # A_hat: perturbed adjacency matrix
  # ClustersLength: length of the k clusters
  # graphs (only if print.plot=TRUE)

  p1 <- p$p_inside
  p2 <- p$p_outside

  if (k>n){
    stop("The number of clusters is larger than the number of nodes. Please choose a smaller number of clusters.")
  }

  MaxLength <- n - 2*(k-1)
  MinLength <- 2
  if (MinLength!=MaxLength){
    Length <- sample(MinLength:MaxLength)
  } else {
    Length <- 2
  }

  if (is.null(ClustersLength)){
    ClustersLength <- c()
    Length_tmp <- Length
    for (i in (1:(k-1))){
      Length_tmp <- Length_tmp[Length_tmp<(n-sum(ClustersLength)+2*(i-k)+1)]
      if (length(Length_tmp)>0){
        clustersLength <- Length_tmp[1]
        ClustersLength <- c(ClustersLength,clustersLength)
        Length_tmp <- Length_tmp[-1]
      } else {
        ClustersLength <- c(ClustersLength,2)
      }
    }
    ClustersLength <- c(ClustersLength,(n-sum(ClustersLength)))
    ClustersLength <- sort(ClustersLength)
  }

  print(paste0(c("There are",k,"clusters of size ",ClustersLength,"."), collapse=" "))

  A <- matrix(1,ncol=ClustersLength[1],nrow=ClustersLength[1])
  for (i in 2:k){
    A <- bdiag(A,matrix(1,ClustersLength[i],ClustersLength[i]));
  }
  A <- matrix(A,ncol(A),nrow(A))

  A <- A-diag(nrow(A)) # adjacency matrix

  # perturbed versions of the graph
  A_perturbed <- matrix(0,n,n)
  edges_inside = edges_outside = 0
  bern_inside = bern_outside <- c()
  for (i in (1:k)){
    bern_inside <- rbinom(ClustersLength[i]*(ClustersLength[i]-1)/2,size=1,prob=1-p1)
    if (i>1){
      bern_outside <- rbinom(sum(ClustersLength[1:(i-1)])*ClustersLength[i],size=1,prob=p2)
    }

    B <- matrix(0,ClustersLength[i],ClustersLength[i])
    B[upper.tri(B, diag = FALSE)] <- bern_inside

    if (i>1){
      A_perturbed[(sum(ClustersLength[1:(i-1)])+1):sum(ClustersLength[1:i]),(sum(ClustersLength[1:(i-1)])+1):sum(ClustersLength[1:i])] <- B
      A_perturbed[(1:sum(ClustersLength[1:(i-1)])),(sum(ClustersLength[1:(i-1)])+1):(sum(ClustersLength[1:i]))] <- bern_outside
    } else {
      A_perturbed[(1:ClustersLength[1]),(1:ClustersLength[1])] <- B
    }
    edges_inside <- edges_inside + length(which(bern_inside==0))
    edges_outside <- edges_outside + length(which(bern_outside==1))
  }
  A_perturbed <- A_perturbed+t(A_perturbed)

  print(paste0("On the ",length(which(A==1))/2," existing edges, ",edges_inside," were removed and ",edges_outside," were added."))

  par(mfrow=c(1,2))
  graph_hat <- graph_from_adjacency_matrix(A_perturbed,mode="undirected")
  graph <- graph_from_adjacency_matrix(A,mode="undirected")
  plot(graph,main="Real graph")
  plot(graph_hat,main="Perturbed graph")

  data <- list(A=A, A_hat = A_perturbed, ClustersLength=ClustersLength)
}
