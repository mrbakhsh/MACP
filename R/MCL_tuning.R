  #' MCL_tuning
  #' @title MCL Hyperparameters Tuning
  #' @param hc_ppi Interactions data containing id1-id2-weight triplets.
  #' @param predcpx A data.frame containing predicted modules resulted from
  #' \code{\link[DeepiCE]{get_clusters}}.
  #' @param refcpx A list containing reference complexes
  #' (i.e., corum complexes).
  #' @param inflation  A vector of integer, representing
  #' MCL inflation parameter
  #' @param csize  An integer, the minimum size of the predicted complexes.
  #' Defaults to 2.
  #' @return  A data.frame containing clustering performance across different
  #' inflation values.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function optimize the choice of MCL algorithm
  #' parameter (inflation) by comparing clustering-derived partitions for each
  #' paramter values to known labels (i.e., CORUM complexes) and
  #' assess the similarity between them using quality measures including
  #' overlap score, sensitivity (Sn),
  #' clustering-wise positive predictive value (PPV), geometric accuracy (Acc),
  #' and maximum matching raio (MMR). It is recommended to first reduce
  #' redundancy in the known reference complexes
  #' via \code{\link{EliminateCpxRedundance}}, then performs parameter tuning.
  #' @export
  #' @examples
  #' # open high-confidence network
  #' hc_ppi <-
  #' read.delim(
  #' system.file("extdata/ppi_input_ClusterONE.txt", package = "MACP"),
  #' header = FALSE)
  #' # predict complexes by ClusterONE
  #'predcpx <-
  #'  get_clusters(csize = 3, d = 0.3, p = 2,
  #'               max_overlap = 0.8,
  #'              tpath = file.path(system.file("extdata", package = "MACP")))
  #'# load the reference complexes
  #'data("refcpx")
  #'# Perform MCL tuning
  #'MCL_tuning_result <-
  #'  MCL_tuning(hc_ppi,
  #'            predcpx,
  #'            refcpx,
  #'           inflation = c(6,8,9))




  MCL_tuning <-
    function(hc_ppi,
             predcpx,
             refcpx,
             inflation = c(6,8,9),
             csize = 2){

      . <- NULL

      if(!is.data.frame(hc_ppi)){
        hc_ppi <- as.data.frame(hc_ppi)
      }

      if(!is.list(refcpx)){
        stop("Reference complexes must be list")
      }
      refcpx <- unname(refcpx)

      colInput <-
        c("ClustID", "Members")
      if(!all(colInput %in% colnames(predcpx))){
        missingCol <-
          setdiff(colInput,
                  colnames(predcpx)[match(colInput, colnames(predcpx))])
        stop("Input data missing: ", paste(missingCol, collapse = ", "))
      }

      if(!is.data.frame(predcpx)){
        predcpx <- as.data.frame(predcpx)
      }
      if(!is.data.frame(hc_ppi)){
        hc_ppi <- as.data.frame(hc_ppi)
      }



      cpx <- list()
      for (i in seq_along(inflation)) {
        cpx[[i]] <-
          MCL_clustering(hc_ppi,predcpx, inflation = inflation[[i]],
                         csize = csize)
      }

      #convert them to the list
      p <- list()
      for(i in seq_along(cpx)) {
        p[[i]] <- str_split(cpx[[i]]$Members, " ")
        names(p[[i]]) <- cpx[[i]]$ClustID
        p[[i]] <-  p[[i]][lapply( p[[i]], length)>=3]
      }
      names(p) <- as.list(inflation)


      # performance assessment
      SCORED_f <- list()
      for(i in seq_along(p)) {
        SCORED_f[[i]] <-  Clust_Valid(p[[i]], refcpx)
      }


      scored_df <-
        data.frame(matrix(unlist(SCORED_f),
                          nrow= length(SCORED_f), byrow=TRUE),
                   stringsAsFactors=FALSE)
      colnames(scored_df) <-
        c("PPV", "Sn", "Acc", "Overlap", "MMR")


      df_out <-
        cbind(as.data.frame(inflation), scored_df) %>%
        mutate(compScore = rowSums(.[,4:6]))
      df_out <- na.omit(df_out)

      return(df_out)
    }


