#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet glmnet
#' @importFrom cvTools cvFolds
#' @title Solve the internal minimization problem
#'
#' @description This internal function of the l1-spectral clustering algorithm solves the l1-minimization problem and recover the community indicators of the clusters.
#' @seealso \code{\link{l1_spectralclustering}}, \code{\link{l1_spectral}}, \code{\link{l1spectral}}.
#' @param U The eigenvector matrix.
#' @param n The number of nodes of the connected component to cluster.
#' @param elements The representative elements of the connected component to cluster.
#' @param iteration The cluster we aim at recovering.
#' @param pen The penalty (to be chosen among "lasso" and "thresholdedLS").
#' @param k The number of clusters.
#'
#' @return \code{v} The community indicator of cluster \code{iteration}.
#' @export
#' @keywords internal
#' @author Camille Champion, Magali Champion
#'
#' @examples
#'  ###################################
#'  # Solving the minimization problem
#'  ###################################
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
#'  n <- ncol(A_tmp)
#'  k <- Clusters$nbr_clusters$Component1 # number of clusters to create
#'  Elements_tmp <- Elements$indices$Component1 # the elements of the first component
#'
#'  # 3rd: perform svd
#'  svd <- eigen(A_tmp)
#'  eigenvalues <- sort(svd$values,index.return=TRUE)
#'  eigenvectors <- svd$vectors[,eigenvalues$ix]
#'
#'  # 4th: solve the minimization problem
#'  i <- 1 # the cluster we aim at recovering
#'  U <- t(eigenvectors[,1:(n-k+i-1)])
#'  v <- PenOpt(U, n, elements = Elements_tmp, iteration = i, pen = "lasso", k) # for lasso
#'  v <- PenOpt(U, n, elements = Elements_tmp, iteration = i, pen = "thresholdedLS", k)

PenOpt <- function(U, n, elements, iteration, pen, k){
  # solve the l1-min problem using Lasso or thresholded least-squares penalization
  w <- U[,elements[iteration]]
  W <- matrix(U[,-elements[iteration]],ncol=(ncol(U)-1),nrow=nrow(U))

  if (sum(w)==0 || length(w)==1){
    print("There is only one node in this cluster.")
    # w is constant - no more nodes in the cluster
    v <- rep(0,n)
    v[elements] <- 1
  } else if (length(which(w!=0))==1){
    sol <- glmnet(W,-w,lambda=0,lower.limits=0)
    sol <- sol$beta
    solution <- as.matrix(sol)

    solution[solution<0.5] <- 0

    if(elements[iteration]==1){
      v <- c( 1 ,solution[elements[iteration]:length(solution)])
    } else{
      v <- c(solution[1:(elements[iteration]-1)], 1 ,solution[elements[iteration]:length(solution)])
    }

  } else {
    if (pen == "lasso"){
      lassosol <- cv.glmnet(W,-w)

      cvtop<- min(lassosol$cvm)+5*(max(lassosol$cvm)-min(lassosol$cvm))/100
      plot(lassosol)
      abline(h=cvtop)

      error <- min(lassosol$cvm[lassosol$cvm>cvtop])
      lambda_opt <- lassosol$lambda[which(lassosol$cvm==error)]

      sol <- glmnet(W,-w,lambda=lambda_opt)
      sol <- sol$beta
      solution <- as.matrix(sol)

      I <- which(solution!=0)
      if (length(I)>1){
        sol2 <- glmnet(W[,I],-w,lambda=0)
        sol2 <- sol2$beta
        sol2 <- as.matrix(sol2)
        solution[I] <- sol2
      }
    } else {
      # thresholded least-squares penalty
      # cross validation test
      if (nrow(W)>=5){
        K <- 5
      } else {
        K <- nrow(W)
      }
      groups <- cvFolds(nrow(W),K=K)
      T <- seq(from=0, to=0.5,0.001)
      erreur_cv <- c()
      for (g in (1:K)){
        Erreur <- c()
        sol <- glmnet(W[which(groups$which!=g),],-w[which(groups$which!=g)],lambda=0)
        sol <- sol$beta
        sol <- as.matrix(sol)

        for (t in (1:length(T))){
          sol2 <- sol
          sol2[which(sol<T[t])] <- 0
          if (length(which(sol2>0))>1){
            solution <- glmnet(W[which(groups$which!=g),which(sol2>0)],-w[which(groups$which!=g)],lambda=0)
            solution <- solution$beta
          } else if (length(which(sol2>0))==1){
            solution <- lm(-w[which(groups$which!=g)]~W[which(groups$which!=g),which(sol2>0)]-1)
            solution <- solution$coefficients
          } else {
            solution <- 0
          }
          sol2[which(sol2>0)] <- solution

          erreur <- sum((-w[which(groups$which==g)]-W[which(groups$which==g),]%*%sol2)^2)
          Erreur <- c(Erreur,erreur)
        }

        erreur_cv <- cbind(erreur_cv,Erreur)
      }
      Erreur_cv <- rowMeans(erreur_cv)

      plot(T,Erreur_cv,type="l")
      lambda_opt <- T[which.min(Erreur_cv)]
      sol2 <- glmnet(W,-w,lambda=0)
      sol2 <- sol2$beta
      sol2 <- as.matrix(sol2)
      sol2[which(sol2<lambda_opt)] <- 0
      solution <- sol2
    }
    if(elements[iteration]==1){
      v <- c( 1 ,solution[elements[iteration]:length(solution)])
    } else if (elements[iteration]==n){
      v <- c(solution[1:(elements[iteration]-1)], 1)
    } else {
      v <- c(solution[1:(elements[iteration]-1)], 1 ,solution[elements[iteration]:length(solution)])
    }
  }
  v[v<0] <- 0
  return(v)
}
