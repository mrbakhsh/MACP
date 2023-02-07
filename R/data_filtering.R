  #' data_filtering
  #' @title Data Filtering
  #' @param x A data matrix  object with rows
  #' including proteins and fractions along the columns.
  #' @return Filtered matrix.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom stats sd
  #' @description This function removes proteins for which peptide only
  #' detected in one fraction (i.e., "one-hit-wonders") across the
  #' co-elution table, common contaminants (e.g., keratins) only for mouse
  #' and human organisms and frequent flyers.
  #' @export
  #' @importFrom utils read.csv
  #' @examples
  #' # Load the co-elution data
  #' data("exampleData")
  #' # Perform raw data pre-processing
  #' datOut <- data_filtering(exampleData)


  data_filtering <- function(x) {


    if (!is.matrix(x)) {
      x <- as.matrix(x)
    }

    if(is.character(x) == TRUE){
      stop("matrix must include numerical variables")
    }

    if(is.null(row.names(x))){
      stop("Please specify the row.names")
    }

    # Remove 'one-hit-wonders'
    x <-
      x[which(rowSums(x != 0) >= 2), ]

    # Remove keratin
    keratin <- read.csv(
      system.file("extdata/kerain_protein.csv", package = "MACP"),
      header = TRUE
    )

    x <- subset(x, !row.names(x) %in% keratin)

    # Remove frequent-flyers
    n <- round(0.8 * ncol(x))
    x <-
      x[which(rowSums(x != 0) < n), ]

    return(x)

  }


