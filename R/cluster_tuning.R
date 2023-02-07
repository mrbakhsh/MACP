
  .out_list <- function(x){
    s <- as.data.frame(x)
    s[,1] <- str_replace_all(s[,1], "\t", " ")
    s <- str_split(s[, 1], " ")
    return(s)

  }


  #' cluster_tuning
  #' @title ClusterONE Hyperparameters Tuning
  #' @param refcpx A list containing reference complexes
  #' (i.e., corum complexes).
  #' @param csize  An integer, the minimum size of the predicted complexes.
  #' Defualts to 2.
  #' @param d A vector of number, density of predicted complexes.
  #' @param p A vector of integer, penalty value for the inclusion of each node.
  #' @param max_overlap A vector of number, specifies the maximum allowed
  #' overlap between two clusters.
  #' @param tpath A character string indicating the path to the project
  #' directory that contains the interaction data. Interaction data must be
  #' stored as .txt file and containing id1-id2-weight triplets. Defaults
  #' to MACP/inst/extdata directory.
  #' @return  A data.frame containing clustering performance across different
  #' combination of parameters.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @references Nepusz, T., Yu, H., and Paccanaro, A. (2012a).
  #' Detecting overlapping protein complexes in protein-protein interaction
  #' networks. Nat. Methods 9, 471.
  #' @importFrom stringr str_replace_all
  #' @importFrom stringr str_split
  #' @description This function optimizes the choice of ClusterONE algorithm
  #' parameters such as density, node penalty, and overlap score by comparing
  #' clustering-derived partitions for each combination of parameters to known
  #' labels (i.e., CORUM complexes) and assess the similarity between them
  #' using quality measures including overlap score, sensitivity (Sn),
  #' clustering-wise positive predictive value (PPV), geometric accuracy (Acc),
  #' and maximum matching raio (MMR).It is recommended to first reduce
  #' redundancy in the known reference complexes via
  #' \code{\link{EliminateCpxRedundance}},then performs parameter tuning.
  #' @export
  #' @importFrom utils write.table
  #' @importFrom tidyr unite
  #' @examples
  #' # first load the reference complex
  #' data("refcpx")
  #' Clust_tuning_result <-
  #' cluster_tuning(refcpx, csize = 3,
  #' d = c(0.3,0.4),
  #' p = c(2, 2.5),
  #' max_overlap = c(0.6,0.7),
  #' tpath = file.path(system.file("extdata", package = "MACP")))






  cluster_tuning<-
    function(refcpx,
             csize = 2,
             d = c(0.3, 0.5),
             p = c(2),
             max_overlap = c(0.5,0.6),
             tpath =
               file.path(system.file("extdata", package = "MACP"))){
      . <- NULL
      name <- NULL


      if(!is.list(refcpx)){
        stop("Reference complexes must be list")
      }
      refcpx <- unname(refcpx)


      # set directory to java file
      fpath <- file.path(system.file("java", package = "MACP"))
      tpath <- tpath


      # create paramter combinations
      max_overlap = max_overlap
      d = d
      p = p

      df <- expand.grid(max_overlap, d, p)
      colnames(df) <- c("max_overlap", "d", "p")
      df$max_overlap <- paste0("--max-overlap ", df$max_overlap)
      df$d <- paste0("-d ", df$d)
      df$p <- paste0("--penalty ", df$p)
      df <- unite(df, name , c(max_overlap, d, p),
                       sep = "~", remove = FALSE)




      txt_EX <-
        paste("java -jar", paste0(fpath,"/","cluster_one-1.0.jar"),
              "--seed-method",
              "nodes",
              df$max_overlap,
              df$d, df$p, "-s", csize, "-F plain",
              paste0(tpath, "/","ppi_input_ClusterONE.txt"))

      javaOutput <-
        lapply(txt_EX, function(x) system(x, intern = TRUE,
                                          ignore.stderr = TRUE))


      listcpx  <- lapply(javaOutput, function(x) .out_list(x))
      names(listcpx) <- df$name
      listcpx <- listcpx[lapply(listcpx, length)>=3]


      SCORED_f <- list()
      for(i in seq_along(listcpx)) {
        SCORED_f[[i]] <-  Clust_Valid(listcpx[[i]], refcpx)
      }
      names(SCORED_f) <- names(listcpx)


      scored_df <-
        data.frame(matrix(unlist(SCORED_f),
                          nrow= length(SCORED_f), byrow=TRUE),
                   stringsAsFactors=FALSE)
      colnames(scored_df) <-
        c("PPV", "Sn", "Acc", "Overlap", "MMR")

      names <-
        data.frame(matrix(unlist(names(SCORED_f)),
                          nrow= length(SCORED_f), byrow=TRUE),
                   stringsAsFactors=FALSE)
      colnames(names)[1] <-
        c("tune_names")

      scored_df <-
        data.frame(matrix(unlist(SCORED_f),
                          nrow= length(SCORED_f), byrow=TRUE),
                   stringsAsFactors=FALSE)
      colnames(scored_df) <-
        c("PPV", "Sn", "Acc", "Overlap", "MMR")

      df_out <-
        cbind(names, scored_df) %>%
        mutate(compScore = rowSums(.[,4:6]))
      df_out <- na.omit(df_out)

      return(df_out)
    }


