#########################################################################################
#   CENTRE FOR BIOINFORMATICS, BIOMARKER-DISCOVERY & INFORMATION-BASED MEDICINE
#   THE UNIVERSITY UNIVERSITY OF NEWCASTLE
#   University Drive, Callaghan, NSW, 2308, AUSTRALIA
#
#   Contains utility functions frequently used in CIBM
#
#   Created on: 2014/02/27
#   Last modified on: 2017/10/19
#   author: Carlos Riveros
#   updated: Renato Vimieiro
#
#   License: MIT <http://opensource.org/licenses/MIT>
#
#   Copyright (c) 2014-2017 Carlos Riveros
#
#  [20140403:RV] Added attribute label to IGRAPH object to make it display
#                the label in yEd.

#' @title       Recursive MSTkNN clustering
#' @name        mstknn
#' @description Runs the MSTkNN iteratively splitting up the clusters based
#'              on the k prescription
#' @details     The distance is used to construct a Minimum Spanning Tree and k-Nearest Neighbors
#'              graphs. All edges in MST not in kNN are removed. The process is repeated
#'              recursively for each resulting connected component, until no more changes in the
#'              graph structure occur. The initial number of neighbors can be provided externally.
#'              By default, it is estimated as
#'              \code{min{floor(log(|V|)), k \ kNN is connected}}.
#'
#'              \bold{TODO:} "Collection" of small clusters (2 - 3 nodes) into larger clusters
#'              based on removed MST edge.
#' @return      A graph object is returned.
#' @export
#' @param d     a matrix or 'dist' object with the distance matrix for the computation.  Should be
#'              symmetrical and positive definite, although this is not checked nor enforced.
#'              Column names, if present, are used to label vertices, otherwise they are imposed as
#'              \code{vNNNN}, where NNNN is the row/column number.
#' @param k     an optional number of neighbors to compute the kNN graph.  The default is
#'              to compute as \code{max{ceil(log(|V|)), k \ kNN is not connected}}.
#'              If given, then no iterative search of a connected kNN graph is performed.
#' @param min.size Minimum cluster size.  Clusters smaller than this size will not be
#'              recursively search for split.
#' @param verbose If messages are desired.
#'
#' @seealso \code{\link[igraph:write_graph]{write_graph}},\code{\link[igraph:components]{components}}
#' @examples
#' require(igraph)
#' a <- matrix(data=runif(20000),nrow=200)
#' colnames(a) <- sapply(1:dim(a)[2], function(x) paste("v",x,sep=""))
#' d <- distance(a)
#' gmstknn <- mstknn(d,verbose=TRUE)
#' components(gmstknn)
#' \dontrun{write_graph(gmstknn,file="test.gml",format="gml")}
#'
#' # Another example with the Ray et al. Alzheimer's data
#' data(alzheimer)
#' d <- distance(alzheimer,method="pearson")
#' gmstknn <- mstknn(d)
#' clusters(gmstknn)
#' \dontrun{write_graph(gmstknn,file="test.gml",format="gml")}

mstknn <- function(d,k=NULL,min.size=10,verbose=FALSE) {
    # Data conversion
    tini <- Sys.time()
    if(class(d) %in% c("dist","data.frame"))    d <- as(d,"matrix")
    stopifnot(inherits(d,'matrix') || inherits(d,'Matrix'))
    # stopifnot(all(d >= 0))
    dmx <- max(d)
    dmn <- min(d)

    tend <- Sys.time()
    if (verbose) {
        cat(sprintf('mstknn: data preprocessing %fs\n',difftime(tend,tini,units="secs")))
        cat(sprintf('Max distance: %f, Min distance: %f\n',dmx,dmn))
    }

    # [CR.20140416] FIXME:  For very large matrices, we should use 'ff'.
    # For now, we just Call the garbage collector with hope.
    gc()
    if(min.size < 4) min.size <- 4
    stopifnot(dim(d)[1] == dim(d)[2])
    nc <- dim(d)[1]
    if(is.null(colnames(d)))
        colnames(d) <- sapply(1:nc, function(x) paste("v",x,sep=""))
    if(is.null(rownames(d)) || any(rownames(d) != colnames(d)))
        rownames(d) <- colnames(d)

    # Adding attribute label to vertices because yEd doesn't recognize 'name'
    .g <- .rmstknn(d,k,min.size=min.size,verbose=verbose)
    V(.g)$label <- V(.g)$name

    return(.g)
}

# Internal function tracking level of recursion
# It asumes all checks on variables have been done on public interface
#
.rmstknn <- function(d,k=NULL,min.size=20,verbose=FALSE,level=0) {
    tini <- Sys.time()
    nc <- dim(d)[1]
    level <- level+1

    if(verbose) {
        cat(sprintf('>> [%d] Start...\n',level))
    }
    # Initial pass
    g <- .mstknn(d=d,k=k,min.size=min.size,verbose=verbose)
    c <- clusters(g)
    if(verbose) {
        cat(sprintf('[%d] Components: %d\n',level,c$no))
        cat('  [',c$csize,' ]\n')
    }
    if(c$no > 1) {
        for(j in 1:c$no) {
            if(c$csize[j] < min.size) {
                if(verbose) {
                    cat(sprintf('[%d] Skip Component %d, size %d (< %d)\n',
                                level,j,c$csize[j],min.size))
                }
                next
            }
            if(verbose) {
                cat(sprintf('[%d] Analysing Component %d, size %d\n',
                                level,j,c$csize[j]))
            }
            s1 <- c$membership == j
            g1 <- .rmstknn(d=d[s1,s1],k=k,min.size=min.size,verbose=verbose,level=level)
            c1 <- clusters(g1)
            if(c1$no > 1) {
                if(verbose) {
                    cat(sprintf("[%d] Split cluster %d into %d\n",level,j,c1$no))
                }
                # Edges no longer included, as character matrix
                ee <- get.edgelist(graph.difference(induced.subgraph(g,s1),g1,byname=TRUE))
                # and then as edge ids on the larger graph... (all this to look for edges by name, puaj)
                ee <- get.edge.ids(g,as.vector(t(ee)))
                if(any(ee == 0))
                    warning("Edge id not found")
                # Remove them
                g <- delete.edges(g,ee)
                # DEBUG
                cqcq <- clusters(g)
                if(verbose) {
                    cat('  [',cqcq$csize,' ]\n')
                }
            }
            else {
                if(verbose) {
                    cat(sprintf("[%d] Cluster %d: no change\n",level,j))
                }
            }
        }
    }
    tend <- Sys.time()
    if(verbose) {
        c <- clusters(g)
        cat(sprintf('[%d] rMSTkNN: %d clusters, %fs\n',level,c$no,difftime(tend,tini,units="secs")))
        summary(g)
        cat('  [',c$csize,' ]\n')
    }
    return(g)
}

# @title MST-kNN clustering
# @name mstknn
# @rdname mstknn
# @author Carlos Riveros
# @description Main step of clustering via the MST-kNN algorithm from (reference).
# @details The distance is used to construct a Minimum Spanning Tree and a k-Nearest Neighbors
# graphs. All edges in MST not in kNN are removed. The process is repeated recursively for each
# resulting connected component, until no more changes in the graph structure occur.
# The initial number of neighbors can be provided externally. By default, it is estimated as
# \code{ceil(|V|)}.
#
#              A graph object is returned.
#
# \bold{TODO:} "Collection" of small clusters (2 - 3 nodes) into larger clusters based on removed
# MST edge.
# @export
# @param d a matrix or 'dist' object with the distance matrix for the computation.  SHould be
#        symmetrical and positive definite, although this is not checked nor enforced.
#        Column names, if present, are used to label vertices.
# @param k an optional number of neighbors to compute the kNN graph.  The default is
#        to compute as \code{max{ceil(log(|V|)), k \ kNN is not connected}}.
#        If given, then no iterative search of a connected kNN graph is performed.
# @param min.size  Size of the smallest cluster to be returned
# @param verbose If messages are desired.
# @seealso \code{\link[stats]{dist}, \link[cluster]{daisy}} and \code{\link["cibm.utils"]{JSD}}
#          for categorical data and other measures
# @return A graph object with the MST-kNN clusters.  Edge weights are distances.
# @import igraph Matrix
# @examples
#
# a <- matrix(data=runif(20000),nrow=200)
# colnames(a) <- sapply(1:dim(a)[2], function(x) paste("v",x,sep=""))
# d <- distance(a)
# gmstknn <- mstknn(d)
# \dontrun{
# plot(gmstknn,layout=layout.reingold.tilford(gmst))
# write_graph(gmstknn,file="random.cluster.gml",format="gml")
# }
.mstknn <- function(d,k=NULL,min.size=10,verbose=FALSE) {
    tini <- Sys.time()
    nc <- dim(d)[1]
    iter <- is.null(k)
    if(iter)
        k <- floor(log(nc))

    # Construct Minimum Spanning Tree
    ai <- Sys.time()
    ogmst <- minimum.spanning.tree(
        graph.adjacency(d,mode='undirected',weighted=TRUE,diag=FALSE),
        algorithm='prim'
    )
    ae <- Sys.time()
    if(verbose) {
        cat(sprintf('  MST: time=%fs\n',difftime(ae,ai,units="secs")))
        summary(ogmst)
    }
    gc()

    # Adjacency matrix for MST
    mmst <- get.adjacency(ogmst,attr='weight')

    r <- .kNNadjacency(d,k,iter,verbose)

    mmst <- mmst * r$adjacency
    ngmst <- graph.adjacency(mmst,mode='undirected',weighted=TRUE,diag=FALSE)
    ngmst <- .joinSmall(ngmst,ogmst,min.small=min.size,verbose=verbose)
    tend <- Sys.time()
    if(verbose) {
        cat(sprintf('  MSTkNN: time=%fs\n',difftime(tend,tini,units="secs")))
        summary(ngmst)
    }
    return(ngmst)
}

# title         .kNNadjacency
# description   Connectivity and adjacency matrix of kNN graph
# details       Compute the adjacency matrix and connectivity of the kNN graph
#               determined by the distance matrix d, starting at the k neighbour nbr.
#               iniK.  The default behavior will iterate down the iniK until the graph
#               becomes disconnected.  Iterate == FALSE means a single step ;-)
#
# param d       The distance matrix
# param iniK    The initial (maximum) value of neighbors
# param iterate Wether to iterate looking for the minimum k such that the kNN is connected
# param grow.k  If no connected graph with up to iniK neigbohrs can be constructed, keep increasing
#               k for at most this value.
#
# return        a list with the following elements:
# section Value:
# \describe{
#  \item{\code{connected}}{A logical indicating if the graph is connected}
#  \item{\code{k}}{The final k value}
#  \item{\code{adjacency}}{A symmetric sparse Matrix with the kNN adjacency}
# }
#
# The computation is not particularly optimised
#' @import Matrix
.kNNadjacency <- function(d, iniK, iterate=TRUE, verbose=FALSE, grow.k=3) {
    nc = dim(d)[1]
    am <- Matrix(0,nrow=nc,ncol=nc,sparse=TRUE)
    colnames(am) <- colnames(d)
    rownames(am) <- colnames(d) # Or it will not be symmetric !

    # om <- matrix(0,nrow=nc,ncol=iniK)   # matrix of indexes
    onu <- matrix(1:nc,nrow=nc,ncol=2)  # To index the am matrix
    onl <- onu

    # TODO: optimise
    # Use row order in matrix to resolve ties in distance
    myf <- function(u,ik=(iniK+grow.k)) order(u,seq(u),na.last=T)[2:(ik+1)]
    aini <- Sys.time()
    ai <- aini
    om <- t(apply(d,1,myf))  # matrix of indexes
#     for(j in 1:nc) {
#         o <- order(d[j,],na.last=T)[1:(iniK+1)]    # First element should be == j
#         if(o[1] != j) message('Another 0 distance: ',j,o[1])
#         om[j,] <- o[2:(iniK+1)]
#     }
    ae <- Sys.time()
    if(verbose) {
        cat(sprintf('  NN: nodes=%d, initial k=%d, time=%fs\n',
                    nc,iniK,difftime(ae,ai,units="secs")))
    }

    isc <- FALSE
    if(!iterate) {
        for(j in 1:iniK) {
            onu[,2] <- om[,j]
            onl[,1] <- om[,j]
            am[onu] <- 1
            am[onl] <- 1
        }
        am <- forceSymmetric(am,uplo='U')

        ai <- Sys.time()
        g1 <- graph.adjacency(adjmatrix=am,mode='undirected',weighted=TRUE,diag=FALSE)
        isc <- is.connected(g1)
        ae <- Sys.time()
        # if(verbose) message('  Check k=',iniK,' ',isc,' ',difftime(ae,ai,,units="secs"))
        result <- list(connected=isc, k=iniK, adjacency=am, knnIdx=om)
    }
    else {
        thisK <- 0
        while(thisK < (iniK+grow.k) && !isc) {
            thisK <- thisK + 1
            ai <- Sys.time()
            onu[,2] <- om[,thisK]
            onl[,1] <- om[,thisK]
            am[onu] <- 1
            am[onl] <- 1
            am <- forceSymmetric(am,uplo='U')
            ae <- Sys.time()
            ai <- Sys.time()
            g1 <- graph.adjacency(adjmatrix=am,mode='undirected',weighted=T,diag=F)
            isc <- is.connected(g1)
            ae <- Sys.time()
            # if(verbose) message('  Check k=',thisK,' ',isc,' ',difftime(ae,ai,units="secs"))
        }
        result <- list(connected=isc, k=thisK, adjacency=am, knnIdx=om)
    }
    if(verbose) {
        ae <- Sys.time()
        cat(sprintf('  NN: final k=%d, conn=%d, time=%fs\n',
                    result$k, result$connected, difftime(ae,aini,units="secs")))
    }
    return(result)
}

# title         joinSmall
# description   Joins small components to larger clusters, based on distance
# details       For components of less than the given size, they will be reconnected to
#               the component closer to any of their vertices via an edge that existed before the split.
#               Joining is applied to smaller components first, and components recalculated.
# param newg    Graph to be searched for small components.  Vertex names and number should be same as
#               vertex names in `oldg`.
# param oldg    Original graph from where `newg` was derived. Distance is obtained from the graph 'weight'
#               edge attribute.
# param min.small  Components of this size and smaller are joined to larger components.
# param verbose
#
# returns       A joined graph based on \code{g}.  New edges are added with a 'weight' attribute
#               containing the distance
.joinSmall <- function(newg, oldg, min.small=4, verbose=T) {
    stopifnot(vcount(newg) == vcount(oldg))   # Until we are sure
    c <- components(newg)
    if(vcount(newg) <= min.small || c$no == 1)
        return(newg)
    if(verbose) {
        cat('Before reduction: V size',vcount(newg),',',c$no,'clusters [',c$csize,']\n')
    }
    while(c$csize[which.min(c$csize)] < min.small) {
        sidx <- which.min(c$csize)
        nvi <- which(c$membership == sidx)
        # Only works if node ids are same !
        ee <- E(oldg)[ nvi %--% V(oldg)[-nvi] ]
        ilnk <- which.min(ee$weight)
        if(verbose) {
            cat('Joining small cluster',sidx,', size',c$csize[sidx],'by link:\n')
            print(ee[ilnk])
        }
        # Convoluted way to get names in old graph to apply to new graph
        newg <- newg + edge(as_edgelist(oldg)[as.numeric(ee[ilnk]),],weight=ee[ilnk]$weight)
        c <- components(newg)
    }
    return(newg)
}

