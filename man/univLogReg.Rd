\name{univLogReg}
\alias{univLogReg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Univariate logistic regression with simplex optimizer
}
\description{
This function implements univariate logistic regression, with a single continuous or ordinal predictor. It is intended for usage as the weak classifier for the adaboost function.
}
\usage{
univLogReg(weights, data, column)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weights}{
  A vector of weights. Each weight should correspond to an observation in the data.frame (the number of rows). Following the adaboost algorithm, misclassified observations are given more weight on following classifications.
  }
  \item{data}{
  A matrix containing the data. The last column should be the outcome variable.
  }
  \item{column}{
  An integer indicating which column of the matrix contains the predictor of interest. The index for this function starts at 0, not 1, since it is a C++ function. Again, note that this function only accepts a single predictor at a time.
  }
}
\details{
This function uses the simplex optimizer shown in the lecture slides. The optimizer is used to maximize the weighted MLE.
}
\value{
Returns a vector<bool> with two entries, containing the coefficients for logistic regression output. The first entry is the intercept, and the second is the coefficient for the predictor.
}
\author{
Dylan Sun and Nina Zhou
}
\note{
This function can only handle one covariate at a time; if you want to use a factor predictor with more than two levels, you will have to manually create dummy columns.
}
\examples{
library(adaboost615)
library(data.table)
library(microbenchmark)
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
train.matrix <- as.matrix(train)

initial.weights <- rep(1/dim(train.matrix)[1], dim(train.matrix)[1])
betas <- univLogReg(initial.weights, train.matrix, 0) # 0 is the first column index
betas <- univLogReg(initial.weights, train.matrix, 1) # 1 is the second column index
microbenchmark(univLogReg(initial.weights, train.matrix, 0))
microbenchmark(univLogReg(initial.weights, train.matrix, 1))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
