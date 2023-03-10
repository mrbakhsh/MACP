% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scaling.R
\name{scaling}
\alias{scaling}
\title{Column and Row-wise Normalization}
\usage{
scaling(data)
}
\arguments{
\item{data}{A data matrix with rows
including proteins and fractions along the columns.}
}
\value{
Scaled data matrix.
}
\description{
This function performs column and row-wise normalization.
}
\details{
scaling
}
\examples{
# Load the co-elution data
data("exampleData")
# Normalize the data
datOut <- scaling(exampleData)
}
\author{
Matineh Rahmatbakhsh, \email{matinerb.94@gmail.com}
}
