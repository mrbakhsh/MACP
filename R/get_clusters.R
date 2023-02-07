  #' get_clusters
  #' @title Predict Complexes
  #' @param csize An integer, the minimum size of the predicted complexes.
  #' Defaults to 2.
  #' @param d A number, density of predicted complexes. Defaults to 0.3.
  #' @param p An integer, penalty value for the inclusion of each node.
  #' Defaults to 2.
  #' @param max_overlap A number, specifies the maximum allowed
  #' overlap between two clusters. Defaults to 0.8.
  #' @param tpath A character string indicating the path to the project
  #' directory that contains the interaction data. Interactions data must be
  #' stored as .txt file and containing id1-id2-weight triplets.
  #' @return  A data.frame containing predicted complexes
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @references Nepusz, T., Yu, H., and Paccanaro, A. (2012a).
  #' Detecting overlapping protein complexes in protein-protein interaction
  #' networks. Nat. Methods 9, 471.
  #' @importFrom tibble rownames_to_column
  #' @importFrom stringr str_replace_all
  #' @description This function partitions high-confidence network to
  #' putative complexes via ClusterONE clustering algorithm to identify
  #' protein complex membership.
  #' @export
  #' @examples
  #' predcpx <-
  #' get_clusters(csize = 3, d = 0.3, p = 2,
  #' max_overlap = 0.8,
  #' tpath = file.path(system.file("extdata", package = "MACP")))

  get_clusters <-
    function(csize = 2,
             d = 0.3,
             p = 2,
             max_overlap = 0.8,
             tpath =
               file.path(system.file("extdata", package = "MACP")))
    {

      # set directory to java file
      fpath <- file.path(system.file("java", package = "MACP"))
      tpath <- tpath

      # Parameter input
      max_overlap = max_overlap
      d = d
      p = p


      txt_EX <-
        paste("java -jar", paste0(fpath,"/","cluster_one-1.0.jar"),
              "--seed-method",
              "nodes",
              "--max-overlap",
              max_overlap,"-d",d, "--penalty", p, "-s", csize,
              paste0(tpath,"/","ppi_input_ClusterONE.txt"))

      javaOutput <-
        system(txt_EX, intern = TRUE, ignore.stderr = TRUE)

      df_clust <- as.data.frame(javaOutput)
      df_clust[,1] <-
        str_replace_all(df_clust[,1], "\t", " ")
      df_clust <-
        rownames_to_column(df_clust, "ClustID")
      colnames(df_clust)[2] <- "Members"

      return(df_clust)
    }

