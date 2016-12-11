\name{errorRate}
\alias{errorRate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Compares two columns and returns the error rate.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is a simple function used to calculate error rate. It compare two columns and returns the percentage of entries that are not equal in the both of them.
}
\usage{
predictAda615(test.data, beta.weights)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{test.data}{
  A matrix containing the test data. The matrix should not include the outcome label column, and should therefore have one less column than the training matrix. The covariates should all otherwise be in the same order as in the training matrix.
  }
  \item{beta.weights}{
  A matrix containing the output from the adaboost function.
  }
}
\details{
Using the coefficients from each iteration of logistic regression, this function generates probabilities with the equation prob = 1/(1 + exp(-1*(train_results(j,0) + train_results(j,1)*test(i,j)))). It then uses 0.5 as the breakpoint for assigning labels.
}
\value{
Returns a vector containing the predict outcome labels. The vector will be the same length as the number of rows in the input matrix.
}
\author{
Dylan Sun and Nina Zhou
}
\note{
The support function errorRate can be then used to calculate test error, or test error can be calculated manually.
}
\examples{
See ?adaboost for a full example, including the predicted function and errorRate function.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
