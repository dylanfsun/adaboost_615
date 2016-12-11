\name{adaboost}
\alias{adaboost}
\alias{adaboost615}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Adaboost with logistic regression as the weak classifier
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function accepts a matrix of data, and outputs a matrix of coefficients and weights.
}
\usage{
adaboost(data)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
  A matrix containing all relevant predictors and the class labels. The class labels should be binary 0 and 1, and all predictors must be continuous because they will be treated as such. Categorical predictors should be ordinal or dummy variables should be manually created. The last column must necessarily contain the outcome variable. If the matrix does not store the outcome variable in the last column, then it must be rearranged before using the function.
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
The first two columns of the test data-age and sex and the outcome variable survival were subsetted from the Titanic dataset. 
With the two main predictor age and sex, Adaboost function can correctly classify 0.68 of the survival cases. 
Our Adaboost function is robust to some noisy data, we first added 9 columns of random numbers to the data and the prediction accuracy remains the same. 
We can see from the variable weights that age and sex weighted much higher than others. The ability of handle noise is diminishing with increasing amount of noisy data. 
When we add 1000 columns of random numbers the model predicts all test cases to be dead. In this case, the model has an error rate of 24%.    See ?test.data for more details.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# Load in the data
data(test.data)

## Example with two main predictors
# with first 656 rows as the training dataset and last 100 as the test dataset
beta_weights<-adaboost(as.matrix(test.data[1:656,c(1:2,1003)]))
pred<-predictAda615(as.matrix(test.data[657:756,1:2]),beta_weights)
errorRate(pred,test.data[657:756,1003])

## After adding 9 columns of noise
beta_weights<-adaboost(as.matrix(test.data[1:656,c(1:11,1003)]))
pred<-predictAda615(as.matrix(test.data[657:756,1:11]),beta_weights)
errorRate(pred,test.data[657:756,1003])

## Adding 1000 clumns of noise
beta_weights<-adaboost(as.matrix(test.data[1:656,]))
pred<-predictAda615(as.matrix(test.data[657:756,-1003]),beta_weights)
errorRate(pred,test.data[657:756,1003])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
