\name{test.data}
\alias{test.data}
\docType{data}
\title{
test.data
}
\description{
Test data for usage in adaboost package for Biostat 615 project.
}
\usage{data(test.data)}
\format{
  A data frame with 756 observations on the following 1003 variables.
  \describe{
    \item{\code{Age}}{a numeric vector containing the age of the passenger}
    \item{\code{Sex}}{a binary numeric vector, containing 0 for female and 1 for male}
    \item{\samp{1, ..., 1000}}{Numeric vectors containing randomly generated data for testing purposes.}
    \item{\code{Survival}}{a binary numeric vector containing the outcome labels.}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
http://math.ucdenver.edu/RTutorial/titanic.txt
}
\examples{
data(test.data)
str(test.data)

# Using data(test.data) should automatically load in the test.data included with the package.
# Otherwise, the data can be manually downloaded and typed as follows:

library(data.table)
train <- read.csv('http://math.ucdenver.edu/RTutorial/titanic.txt',sep='\t')
train <- data.table(train)
train[, Name := NULL]
train[, PClass := NULL]
train <- train[complete.cases(train)]
train[, Sex := as.character(Sex)]
train[Sex == "female", Sex := "0"]
train[Sex == "male", Sex := "1"]
train[, Sex := as.numeric(Sex)]
train[, Age := as.numeric(Age)]
train[, Survived := as.numeric(Survived)]
}
\keyword{datasets}
