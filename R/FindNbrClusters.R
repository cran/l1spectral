#' @importFrom graphics points
#' @importFrom graphics legend
#' @importFrom graphics abline
#' @title Find the optimal number of clusters
#'
#' @description This internal function of the l1-spectral algorithm finds the optimal number of clusters to build.
#' @param A The adjacency matrix
#' @param structure Output of the function \code{FindStructure()}.
#' @param k True number of clusters (not necessarily needed). If not provided, k is chosen by spectral eigengap.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item{\code{nbr_clusters}}{ Optimal number of clusters by component,}
#' \item{\code{nbr_clusters_total}}{ Optimal total number of clusters.}
#' }
#' @export
#' @seealso \code{\link{l1_spectralclustering}}, \code{\link{l1spectral}}.
#' @author Camille Champion, Magali Champion
#' @examples
#'  #########################################
#'  # Finding the optimal number of clusters
#'  #########################################
#'
#'  # 1st example: non-perturbed graph
#'  Data <- CreateDataSet(k=3, n=20, p=list(p_inside=0,p_outside=0))
#'
#'  Structure <- FindStructure(Data$A_hat)
#'
#'  Clusters <- FindNbrClusters(A = Data$A_hat, structure = Structure, k=3)
#'  # The number of clusters is provided (3): each of the 3 components will be divided into 1 cluster
#'
#'  Clusters <- FindNbrClusters(A = Data$A_hat, structure = Structure, k=5)
#'  # The number of clusters is provided (5) and larger than the number of components (3),
#'  # the spectral eigengap method is used to find the optimal number of clusters of each component.
#'
#'  # 2nd example: perturbed graph
#'  Data <- CreateDataSet(k=3, n=20, p=list(p_inside=0.1,p_outside=0.1))
#'
#'  Structure <- FindStructure(Data$A_hat) # there are less than 3 components
#'
#'  Clusters <- FindNbrClusters(A = Data$A_hat, structure = Structure)
#'  # The number of clusters is optimized using the spectral eigengap method

FindNbrClusters <- function(A, structure, k = NULL){
  # A: the adjacency matrix of the graph
  # structure: already existing connected components
  # k: number of clusters
  # Outputs:
  #   - nbr_clusters: optimal number of clusters by component
  #   - nbr_clusters_total: optimal total number of clusters

  nbr_group <- length(structure$groups) # number of connected components

  Eigen_list <- function(group){
    # spectral decomposition on subgraphs of a graph
    if (length(group)>1){
      D <- diag(rowSums(A[group,group])) # degree matrix
      L <- D-A[group,group]
      svd <- eigen(L)
      eigenvalues <- sort(svd$values)
    } else {
      eigenvalues <- NA
    }
    return(eigenvalues)
  }

  Gap <- function(eigenvalue){
    # compute the spectral gap
    if (length(eigenvalue)==1){
      nbr_cluster <- 1
    } else {
      gap <- c()
      for (i in (2:length(eigenvalue))){
        gap <- c(gap,eigenvalue[i]-eigenvalue[i-1])
      }
      gap_ecart <- c(gap,0)-c(0,gap)
      nbr_cluster <- which(gap_ecart[2:length(gap_ecart)]>0.20)[1]+1
      if (is.na(nbr_cluster)){
        nbr_cluster=1
      }
    }
    return(nbr_cluster)
  }

  # various possible cases
  if (is.null(k)){
    # k is not provided

    eigenvalues <- lapply(X = c(structure$groups,list(all=c(1:ncol(A)))),FUN = Eigen_list)

    par(mfrow=c(1,1))
    plot(eigenvalues$all,main="Eigenvalues of the Laplacian matrix",ylab="Eigenvalues",xlab="",type="b")

    if (nbr_group>1){
      for (i in (1:nbr_group)){
        points(eigenvalues[[i]],col=rainbow(nbr_group)[i],type="b")
      }
      legend("bottomright",legend = c("All nodes",paste0("Connected component ",c(1:nbr_group))),col=c("black",rainbow(nbr_group)[1:nbr_group]),lty=1)
    }

    gaps <- lapply(X = eigenvalues[-length(eigenvalues)],FUN = Gap)
    par(xpd=FALSE)
    for (i in (1:length(gaps))){
      color <- c(rainbow(nbr_group),"black")
      abline(v=gaps[[i]],col=color[i],ylim=c(0,10),lty=3)
    }

    if (length(gaps)==1){
      nbr_clusters_total <- gaps[[1]]
      print(paste0("The optimal number of clusters is ",nbr_clusters_total,"."))
    } else {
      nbr_clusters_total <- sum(unlist(gaps))
      print(paste0("The optimal number of clusters is ",nbr_clusters_total,"."))
      print(paste0(c("Here,",nbr_group,"connected components were detected. Each of them should be clustered into",unlist(gaps[1:nbr_group]),"clusters."),collapse=" "))
    }

  } else if (nbr_group==1 && !is.null(k)){
    # easy case: one component and the number of clusters k is provided

    nbr_clusters_total <- k
    gaps <- list(Component1=k)
    print(paste0("The provided number of clusters is ",nbr_clusters_total,"."))

  } else {
    # last case: more than one component and the number of clusters k is provided

    if (k<nbr_group){

      print(paste0("The provided number of clusters ",k," is smaller than the number of components. This is not possible, we thus adjust the number of clusters to the number of components."))
      k <- nbr_group
      nbr_clusters_total <- k
      gaps <- rep(list(1),nbr_group)
      names(gaps) <- paste0("Component",1:nbr_group)
      print(paste0("The total number of clusters is ",k,"."))

    } else if (k==nbr_group){

      k <- nbr_group
      nbr_clusters_total <- k
      gaps <- rep(list(1),nbr_group)
      names(gaps) <- paste0("Component",1:nbr_group)
      print(paste0("The provided number of clusters is ",k,"."))

    } else {

      eigenvalues <- lapply(X = c(structure$groups,list(all=c(1:ncol(A)))),FUN = Eigen_list)
      par(mfrow=c(1,1))
      plot(eigenvalues$all,main="Eigenvalues of the Laplacian matrix",ylab="Eigenvalues",xlab="",type="b")
      for (i in (1:nbr_group)){
        points(eigenvalues[[i]],col=rainbow(nbr_group)[i],type="b")
      }
      legend("bottomright",legend = c("All nodes",paste0("Connected component ",c(1:nbr_group))),col=c("black",rainbow(nbr_group)[1:nbr_group]),lty=1)

      gaps <- lapply(X = eigenvalues[-length(eigenvalues)],FUN = Gap)
      par(xpd=FALSE)
      for (i in (1:length(gaps))){
        color <- c(rainbow(nbr_group),"black")
        abline(v=gaps[[i]],col=color[i],ylim=c(0,10),lty=3)
      }

      if (sum(unlist(gaps)[1:(length(gaps))]) != k){
        # there is a problem
        message(paste0("The optimal number of clusters ",sum(unlist(gaps)[1:(length(gaps))])," does not coincide to the provided number of clusters."))
        gaps[[(length(gaps))]] <- k-sum(unlist(gaps)[1:(length(gaps)-1)])
      }
      nbr_clusters_total <- sum(unlist(gaps)[1:(length(gaps))])
      print(paste0("The provided number of clusters is ",k,"."))
      print(paste0(c("Here,",nbr_group,"connected components were detected. Each of them should be clustered into",unlist(gaps[1:nbr_group]),"clusters."),collapse=" "))
    }

  }

  return(list(nbr_clusters=gaps,nbr_clusters_total=nbr_clusters_total))
}
