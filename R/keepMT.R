  #' keepMT
  #' @title Keep Mitochondrial (mt) Proteins
  #' @param x A data matrix  object with rows
  #' including proteins and fractions along the columns.
  #' @return Matrix containing mt proteins.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @description This function removes all the non-mitochondrial proteins by
  #' mapping the co-eluted proteins from chromatography fractions to
  #' MitoCarta database.
  #' Note that this function is only applicable to mouse or human organisms.
  #' @export
  #' @examples
  #' # Load the co-elution data
  #' data("exampleData")
  #' # Removes non-mitochondtial proteins
  #' datOut <- keepMT(exampleData)


   keepMT <- function(x) {


    if (!is.matrix(x)) {
      x <- as.matrix(x)
    }

    if(is.character(x) == TRUE){
      stop("matrix must include numerical variables")
    }

    if(is.null(row.names(x))){
      stop("Please specify the row.names")
    }

    # Keep mt
    mt_proteins <- read.csv(
      system.file("extdata/MitoCarta_protein.csv", package = "MACP"),
      header = TRUE
    )

    x <- subset(x, row.names(x) %in% mt_proteins$p)

    return(x)

  }


