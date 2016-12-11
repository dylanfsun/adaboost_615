\name{predictAda615}
\alias{predictAda615}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Predict outcome labels using adaboost results.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function accepts new test data and the output from adaboost, and uses them to generate predicted outcome labels.
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
## Time complexity for all three cases
microbenchmark(
pred1<-predictAda615(as.matrix(test.data[657:756,1:2]),beta_weights1),
pred2<-predictAda615(as.matrix(test.data[657:756,1:11]),beta_weights2),
pred3<-predictAda615(as.matrix(test.data[657:756,-1003]),beta_weights3)
)

#Results
Unit: microseconds
                                                                       expr       min         lq
   pred1 <- predictAda615(as.matrix(test.data[657:756, 1:2]), beta_weights1)   172.731   190.2745
  pred2 <- predictAda615(as.matrix(test.data[657:756, 1:11]), beta_weights2)   314.540   347.1540
 pred3 <- predictAda615(as.matrix(test.data[657:756, -1003]), beta_weights3) 16892.515 18488.6240
       mean     median         uq       max neval cld
   299.6531   271.0825   358.2105   999.393   100  a 
   543.6753   468.4010   644.3980  1907.540   100  a 
 22545.6946 20089.5670 24562.7120 61370.046   100   b  
  
See ?adaboost for a full example, including the predicted function and errorRate function.
  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
