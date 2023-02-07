    #' scaling
    #' @title Column and Row-wise Normalization
    #' @param data A data matrix with rows
    #' including proteins and fractions along the columns.
    #' @return Scaled data matrix.
    #' @author Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
    #' @description This function performs column and row-wise normalization.
    #' @export
    #' @examples
    #' # Load the co-elution data
    #' data("exampleData")
    #' # Normalize the data
    #' datOut <- scaling(exampleData)





  scaling <- function(data) {

    if (!is.matrix(data)) {
      data <- as.matrix(data)
    }

    if(is.character(data) == TRUE){
      stop("matrix must include numerical variables")
    }

    if(is.null(row.names(data))){
      stop("Please specify the row.names")

    }
    # column wise
    scale_c <- apply(data, 2,
                     function(x) (x)/sum(x,  na.rm = TRUE))
    scale_c[is.na(scale_c)] <- 0
    # row wise
    scale_r <- apply(scale_c, 1,
                     function(x) (x)/sum(x, na.rm = TRUE))
    scale_r <- t(scale_r)
    scale_r[is.na(scale_r)] <- 0


    return(scale_r)

  }
