  #' EliminateCpxRedundance
  #' @title Hierarchical Clustering of Modules
  #' @param rawCpx A list containing protein complexes
  #' @param custom_bg Vector of proteins names to use as a background.
  #' If given, refcpx will be first mapped to the background proteisn, followed
  #' by removing redundancy in the refcpx.
  #' @param sim_method c(euclidean",
  #' "maximum", "manhattan", "canberra", "binary", or "minkowski"); Default
  #' is euclidean
  #' @param linkage c("average", "ward", "single", "complete", "mcquitty",
  #' "median", "centroid"); Default is average.
  #' @param h numeric scalar or vector with heights where the tree should be
  #' cut; Defaults to 0.2
  #' @return List of unique complexes.
  #' @importFrom stats as.dist
  #' @importFrom stats hclust
  #' @importFrom stats cutree
  #' @importFrom utils unstack
  #' @description This function reduces redundancy in the reference complexes
  #' by first computing the overlap of two complexes via Jaccard index,
  #' followed by merging overlapping complexes with user-defined threshold
  #' (here is 0.2).
  #' @author Matineh Rahmatbakhsh
  #' @export
  #' @examples
  #' # predicted interactions
  #' pred_ppi <- read.table(
  #' system.file("extdata/ppi_input_ClusterONE.txt", package = "MACP"),
  #' header = FALSE)
  #' # get all the proteins in the predicted network
  #' custom_bg <- union(pred_ppi$V1, pred_ppi$V2)
  #' # reference complexes
  #' data("refcpx")
  #' # reduce redundancy in reference complexes
  #' filt_cpx <- EliminateCpxRedundance(refcpx,
  #' custom_bg,
  #' sim_method = "euclidean",
  #' linkage="average",
  #' h = 0.2)



  EliminateCpxRedundance<- function(rawCpx,
                                  custom_bg = NULL,
                                  sim_method = "euclidean",
                                  linkage="average",
                                  h = 0.2){

    if(!is.list(rawCpx)){
      stop("Reference complexes must be list")
    }
    rawCpx <- unname(rawCpx)

    if(!is.null(custom_bg)){
      if(!is.vector(custom_bg)){
        stop("Reference complexes must be vector")
      }
      rawCpx <- lapply(rawCpx, function (x) x[x %in% custom_bg])
      rawCpx <- rawCpx[lapply(rawCpx, length)>=3]


    }

    if(is.null(custom_bg)) {

      rawCpx <- rawCpx
    }

    jd <-# Jaccard Index
      lapply(rawCpx, function(x){
        lapply(rawCpx, function(y){
          return((length(intersect(x, y)))/(length(union(x, y))))
        })
      })
    jm <-
      matrix(unlist(jd), byrow = TRUE, ncol = length(jd))


      dist <- dist(jm, method=sim_method)

    jtree <-
      hclust(dist, method = linkage)
    jbranch <-
      cutree(jtree, h = h)
    names(rawCpx) <-
      paste("CID", jbranch, sep = "_")
    refCpx <-
      lapply(unstack(stack(rawCpx)), function(x){
        return(unique(x))
      })

    return(refCpx)
  }
