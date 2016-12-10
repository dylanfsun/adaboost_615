\name{adaboost_615}
\alias{adaboost_615}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Adaboost with logistic regression as the weak classifier
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function reads in a dataframe of data, and creates an object of class "adamodel" that can be used for prediction.
}
\usage{
adaboost_615(data, dependent_index = dim(data)[2])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A data frame containing all relevant predictors and the class labels. The class labels should be binary 0 and 1.
  }
  \item{dependent_index}{
  An optional integer indicating which column of the data dataframe contains the class labels. By default, adaboost_615 will assume that the last column holds the class labels.
  }

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
Returns an object of class "adamodel" that can be used for prediction.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dylan Sun and Nina Zhou
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
