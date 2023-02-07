  #' enrichmentCPX
  #' @title Functional Enrichment Analysis for Predicted Complexes
  #' @param predcpx A data.frame containing predicted complexes resulted from
  #' \code{\link[MACP]{get_clusters}} or \code{\link[MACP]{MCL_clustering}}.
  #' @param threshold Custom p-value threshold for significance.
  #' @param sources A vector of data sources to use.
  #' See \code{\link[gprofiler2]{gost}} for more details.
  #' @param p.corrction.method The algorithm used for multiple testing
  #' correction;defaults to 'bonferroni'.
  #' See \code{\link[gprofiler2]{gost}} for more details.
  #' @param org An organism name;defaults to 'mmusculus'.
  #' See \code{\link[gprofiler2]{gost}} for more details.
  #' @param custom_bg vector of gene names to use as a statistical background.
  #' Defaults to NULL.
  #' @return A data.frame with the enrichment analysis results.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom dplyr arrange
  #' @importFrom stringr str_split
  #' @description This function uses \code{\link[gprofiler2]{gost}} function
  #' in \code{gprofiler2} package to perform functional enrichment analysis
  #' for predicted modules.
  #' @export
  #' @examples
  #' # predict complexesMCL_clustering
  #' predcpx <- get_clusters()
  #' # perform enrichment for KEGG
  #' enrichCPX <-
  #' enrichmentCPX(predcpx,
  #' sources = "KEGG",
  #' org = "mmusculus")


  enrichmentCPX <-
    function(predcpx,
             threshold = 0.05,
             sources = c("GO", "KEGG", "CORUM", "REAC", "CORUM"),
             p.corrction.method = "bonferroni",
             custom_bg = NULL,
             org = "mmusculus") {

      colInput <-
        c("ClustID", "Members")
      if(!all(colInput %in% colnames(predcpx))){
        missingCol <-
          setdiff(colInput,
                  colnames(predcpx)[match(colInput, colnames(predcpx))])
        stop("Input data missing: ", paste(missingCol, collapse = ", "))
      }

      if(!is.data.frame(predcpx)){
        stop("Input data should be data.frame")
      }

      indcpx <-
        str_split(predcpx$Members, "\\s+")
      names(indcpx) <- predcpx$ClustID

      if(!is.null(custom_bg)) {

        annotCov <- lapply(indcpx, function(x) {
          gprofiler2::gost(x,
                           significant = TRUE,
                           exclude_iea = TRUE,
                           evcodes = TRUE,
                           user_threshold = threshold,
                           sources = sources,
                           correction_method = p.corrction.method,
                           organism = org
          )
        })


      }

      if(is.null(custom_bg)) {
      annotCov <- lapply(indcpx, function(x) {
        gprofiler2::gost(x,
                         significant = TRUE,
                         exclude_iea = TRUE,
                         evcodes = TRUE,
                         user_threshold = threshold,
                         sources = sources,
                         correction_method = p.corrction.method,
                         organism = org
        )
      })
      }

      df <-
        lapply(annotCov, function(x) as.data.frame(x[[1]]))
      ans <-
        purrr::map_df(df, ~ as.data.frame(.x), .id = "id")

      return(ans)
    }
