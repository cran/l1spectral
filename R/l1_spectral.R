#' @title Run the l1-spectral clustering algorithm on one component
#'
#' @description This function runs the l1-spectral clustering algorithm on one component only.
#' @param A The adjacency matrix of the graph to cluster.
#' @param k The number of clusters.
#' @param elements The representative elements of the connected component to cluster.
#' @param pen The penalty (to be chosen among "lasso" and "thresholdedLS").
#' @param stab TRUE/FALSE indicating whether the representative elements should be stabilized (TRUE by default).
#'
#' @return The matrix of community indicators.
#' @export
#' @seealso \code{\link{l1_spectralclustering}}, \code{\link{l1spectral}}.
#' @keywords internal
#' @author Camille Champion, Magali Champion
#' @examples
#'  #########################################################
#'  # Performing the l1-spectral clustering on one component
#'  #########################################################
#'
#'  # 1st: create data
#'  Data <- CreateDataSet(k=3, n=20, p=list(p_inside=0.1,p_outside=0.1))
#'
#'  # 2nd: find the structure, the opt number of clusters and the representative elements
#'  Structure <- FindStructure(Data$A_hat)
#'  Clusters <- FindNbrClusters(A = Data$A_hat, structure = Structure)
#'  Elements <- FindElement(A = Data$A_hat, structure = Structure, clusters = Clusters)
#'
#'  Structure_tmp <- Structure$groups[[1]] # the first component
#'  A_tmp <- Data$A_hat[Structure$groups[[1]],Structure$groups[[1]]]
#'  k <- Clusters$nbr_clusters$Component1 # number of clusters to create
#'  Elements_tmp <- list(score = Elements$score$Component1,
#'                       indices = Elements$indices$Component1)
#'        # the elements of the first component
#'
#'  # 3rd: perform the l1-spectral clustering algorithm
#'  # whithout stabilization
#'  comm <- l1_spectral(A = A_tmp, k = k, elements = Elements_tmp, pen = "lasso", stab=FALSE)
#'
#'  # with stabilization (could take more time)
#'  comm <- l1_spectral(A = A_tmp, k = k, elements = Elements_tmp, pen = "lasso", stab=TRUE)

l1_spectral <- function(A, k, elements, pen, stab = TRUE){
  # A: the adjacency matrix of the graph
  # k: the number of clusters to form (output of the function FindClusters())
  # elements: representative elements of the clusters (output of the function FindElement())
  # pen: the penalty (to be chosen among lasso and threshold)
  # stab: should the indices stabilized? TIME CONSUMING
  # Outputs:
  #   comm: matrix of the components

  if (length(A)==1){
    # only one node in the community
    comm <- 1
  } else {
    # code for running the l1-spectral algorithm for one component
    indices <- elements$indices

    # 1st step: svd on A
    n <- ncol(A)
    svd <- eigen(A)
    eigenvalues <- sort(svd$values,index.return=TRUE)
    eigenvectors <- svd$vectors[,eigenvalues$ix]

    # 2nd step: loop on the number of clusters
    algo <- "stop"
    comm <- c()
    DoubleNodes <- c()

    while (algo == "stop"){
      if (length(DoubleNodes)>0){
        # find other indices
        I <- elements$score[-which(names(elements$score)%in%DoubleNodes)]
        if (length(I)==0){
          print("One cluster disappears.")
          DoubleNodes <- c()
          algo <- "continue"
          break
        } else if (length(I)<k){
          print("One cluster disappears.")
          DoubleNodes <- c()
          comm <- c()
          k <- k-1
          indices <- elements$indices[2:(k+1)]
        } else {
          I <- names(sort(I[1:k]))
          indices <- as.numeric(substring(I, 5))
        }
      }

      if (k>1){
        eigenvectors_tmp <- eigenvectors
        comm <- c()
        for (i in (1:k)){
          # 3rd step: check the indices (only if i>1)
          if (stab==TRUE){
            if (i>1){
              if (length(which(v[indices[-(i-1)]]>0))>0){
                print("Find other community indices.")
                doubleNodes <- paste0("Node",indices[i-1])
                DoubleNodes <- c(DoubleNodes,doubleNodes)
                algo <- "stop"
                break
              } else {
                algo <- "continue"
              }
            }
          }

          # 4th step: Gram-Schmidt (only if i>1)
          if (i>1){
            eigenvectors_tmp <- eigenvectors_tmp[,-(n-k+i-1)]
            eigenvectors_tmp <- cbind(v,eigenvectors_tmp)

            eigenvectors_tmp <- grahm_schmidtCpp(eigenvectors_tmp)
            eigenvectors_tmp <- eigenvectors_tmp$Q
            eigenvectors_tmp <- cbind(eigenvectors_tmp[,2:(n-k+i-1)],eigenvectors_tmp[,1],eigenvectors_tmp[,(n-k+i):n])
          }

          # 5th step: solve the lasso
          U <- t(eigenvectors_tmp[,1:(n-k+i-1)])
          v <- PenOpt(U, n, elements = indices, iteration = i, pen=pen, k)
          print(paste0("Iteration ",i," done."))

          # 6th step: save the community index
          comm <- cbind(comm,v)

          I <- which(rowSums(comm>0)>1)
          if (length(I)>0){
            for (j in (1:length(I))){
              C <- comm[I[j],]
              C[C<max(C)] <- 0
              comm[I[j],] <- C
            }
          }
          if (stab==TRUE){
            if (i==k){
              # check the indices for the last time
              if (length(which(v[indices[-i]]>0))>0){
                print("Find other community indices.")
                doubleNodes <- paste0("Node",indices[-i][which(v[indices[-i]]>0)])
                DoubleNodes <- c(DoubleNodes,doubleNodes)
                algo <- "stop"
                break
              } else {
                algo <- "continue"
              }
            }
          } else {
            algo <- "continue"
          }
        }
      } else {
        comm <- rep(1,n)
        algo <- "continue"
      }
    }
  }
  return(comm)
}
