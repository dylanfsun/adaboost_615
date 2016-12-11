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
  An integer indicating which column of the matrix contains the predictor of interest. Again, note that this function only accepts a single predictor at a time.
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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
