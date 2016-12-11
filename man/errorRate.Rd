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
errorRate(predicted, labels)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{predicted}{
  Column of predicted labels.
  }
  \item{labels}{
  Column of actual labels. Must be the same length as the predicted column.
  }
}
\value{
Returns a double representing the error rate, or proportion of entries that are inequal in the two columns.
}
\author{
Dylan Sun and Nina Zhou
}
\examples{
See ?adaboost for a full example, including the predicted function and errorRate function.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
