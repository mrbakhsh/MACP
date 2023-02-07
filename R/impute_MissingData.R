  #' impute_MissingData
  #' @title Impute missing Values in Elution Profile Matrix
  #' @param x A data matrix with rows
  #' including proteins and fractions along the columns, while some
  #' fractions may contain missing values.
  #' @return Imputed matrix.
  #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
  #' @importFrom zoo na.approx
  #' @description This function imputes missing values in protein elution
  #' profile matrix via average of adjacent rows. This function is
  #' not applicable for missing values present in the first or last column.
  #' @export
  #' @examples
  #' # Load the co-elution data
  #' data("exampleData")
  #' # Replace the values with NAs in the 10th column
  #' exampleData[, 10] <- NA
  #' # Impute missing value
  #' datOut <- impute_MissingData(exampleData)



  impute_MissingData <- function(x) {

    if (!is.matrix(x)) {
      x <- as.matrix(x)
    }

    if(is.character(x) == TRUE){
      stop("matrix must include numerical variables")
    }

    if(is.null(row.names(x))){
      stop("Please specify the row.names")
    }
    clnames <-  colnames(x)

    m_imputed <- t(na.approx(t(x)))
    colnames(m_imputed) <- clnames
    return(m_imputed)
  }
