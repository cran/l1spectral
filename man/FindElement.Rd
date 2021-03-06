% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FindElement.R
\name{FindElement}
\alias{FindElement}
\title{Find the representative elements of the clusters}
\usage{
FindElement(A, structure, clusters, elements = NULL)
}
\arguments{
\item{A}{The adjacency matrix}

\item{structure}{Output of the function \code{FindStructure()}.}

\item{clusters}{Output of the function \code{FindNbrClusters()}.}

\item{elements}{The representative elements of the clusters (not necessary needed). If not provided, chosen using the betweeness centrality score.}
}
\value{
A list with the following elements:
\itemize{
\item{\code{score}}{ The edge betweenness score of all nodes,}
\item{\code{Nodes}}{ Vector of the representative elements.}
}
}
\description{
This internal function of the l1-spectral clustering algorithm finds representative elements of the clusters, that is nodes belonging to the clusters.
}
\examples{
 ######################################################
 # Finding the representative elements of the clusters
 ######################################################

 # 1st: create data (not perturbed graph)
 Data <- CreateDataSet(k=3, n=20, p=list(p_inside=0,p_outside=0))

 # 2nd: find the structure of the graph
 Structure <- FindStructure(Data$A_hat)

 # 3rd: find the optimal number of clusters (here, 3 clusters)
 Clusters <- FindNbrClusters(A = Data$A_hat, structure = Structure, k=3)

 # 4th: find the representative elements of the clusters
 Elements <- FindElement(A = Data$A_hat, structure = Structure, clusters = Clusters)
 # if elements is not provided, the representative elements of each component are chosen
 # by maximizing the edge betweenness score

 Elements <- FindElement(A = Data$A_hat, structure = Structure,
                   clusters = Clusters, elements = c(1,5,12))
}
\seealso{
\code{\link{l1_spectralclustering}}, \code{\link{l1spectral}}.
}
\author{
Camille Champion, Magali Champion
}
