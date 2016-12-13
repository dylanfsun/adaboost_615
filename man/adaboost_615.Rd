\name{adaboost}
\alias{adaboost}
\alias{adaboost615}
\alias{adaboost_615}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Adaboost with logistic regression as the weak classifier
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function accepts a matrix of data, and outputs a matrix of coefficients and weights that result from the adaboost algorithm.
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
See the below Examples section for full R code, without text samples. This section details the results of testing from the Examples code. Alternatively, see the vignette

The first two columns of the test data, Age and Sex, and the outcome variable, Survival, were subsetted from the Titanic dataset. See \code{\link{test.data}} for more details on the original dataset.
\code{data(test.data)}

Using the two main predictors age and sex, the adaboost function can correctly classify 0.68 of the survival cases.

We then tested our adaboost function by adding in several columns of randomly generated noise, starting from 9 extra columns up to 1000 extra columns. Adding 9 columns of random numbers to the data resulted in the same prediction accuracy.

Note that for a fairly small amount of noise columns, the order of the variables does not appear to matter, despite adaboost looping through the dataset in order. When tested on 9 covariates, the model prediction rate remains the same whether the two true predictors, Age and Sex, are placed before or after the randomly generated noise variables.

We can see from the previous example with 9 predictors that age and sex are weighted much more highly than the random noise variables. This ability to handle random noise covariates diminishes as more noisy data is added. When we added 1000 columns of random numbers to the data, the model increased its prediction accuracy, but did so just by predicting all test cases to be dead. In this case, the model has an error rate of 24 percent.

For speed, according to microbenchmark, the adaboost function is faster than logistic function in R when fitting two covariates. See the example for details. The glm function in R failed to coverge when fit on the dataset with 1002 predictors and could not be compared.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
A matrix of coefficients and weights that result from the adaboost algorithm.

Each row of the outcome matrix represents a weak classifier run, e.g. if there are 1002 covariates, adaboost will run 1002 logistic regressions, and there will be 1002 rows in the output matrix. The first column contains the intercept coefficients; the second column contains the predictor coefficients; and the last column contains the weight of each weak classifier.
}
\references{
%% ~put references to the literature/web site here ~
Soft Margins for Adaboost:

http://link.springer.com/article/10.1023/A:1007618119488
}
\author{
Dylan Sun and Nina Zhou
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
See documentation for our supporting functions:

\code{
\link{univLogReg}
\link{predictAda615}
\link{errorRate}
\link{test.data}
}
}
\examples{
# Load in the data
data(test.data) # see ?test.data

## Example with two main predictors
# with first 656 rows as the training dataset and last 100 as the test dataset
beta_weights1<-adaboost(as.matrix(test.data[1:656,c(1:2,1003)]))
pred1<-predictAda615(as.matrix(test.data[657:756,1:2]),beta_weights1)
errorRate(pred1,test.data[657:756,1003])

## After adding 9 columns of noise
beta_weights2<-adaboost(as.matrix(test.data[1:656,c(1:11,1003)]))
pred2<-predictAda615(as.matrix(test.data[657:756,1:11]),beta_weights2)
errorRate(pred2,test.data[657:756,1003])

## Adding 1000 clumns of noise
beta_weights3<-adaboost(as.matrix(test.data[1:656,]))
pred3<-predictAda615(as.matrix(test.data[657:756,-1003]),beta_weights3)
errorRate(pred3,test.data[657:756,1003])

## Time complexity for all three cases
microbenchmark(
beta_weights1<-adaboost(as.matrix(test.data[1:656,c(1:2,1003)])),
beta_weights2<-adaboost(as.matrix(test.data[1:656,c(1:11,1003)])),
beta_weights3<-adaboost(as.matrix(test.data[1:656,]))
)

## Time complexity results
#Unit: milliseconds
#                                                                 expr         min         lq
#  beta_weights1 <- adaboost(as.matrix(test.data[1:656, c(1:2, 1003)]))    2.344701    2.48756
# beta_weights2 <- adaboost(as.matrix(test.data[1:656, c(1:11, 1003)]))   11.339446   12.12715
#              beta_weights3 <- adaboost(as.matrix(test.data[1:656, ])) 1040.429919 1095.18169
#        mean      median          uq         max neval cld
#    2.998985    2.628381    3.268289    8.735227   100  a
#   13.321431   12.748101   13.803661   21.472355   100  a
# 1131.221948 1119.005398 1158.325284 1375.517888   100   b
}

# Compare with logistic regression in R
microbenchmark(
r.log<-glm(Survival~Age+Sex,family = binomial(link = "logit"),data = test.data)
)

#Unit: milliseconds
#
#r.log <- glm(Survival ~ Age + Sex, family = binomial(link = "logit"),
#expr                   min       lq     mean   median       uq      max neval
#data = test.data) 3.418167 3.960986 4.440878 4.241264 4.941011 6.826317   100


% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
