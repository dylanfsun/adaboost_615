#include <Rcpp.h>
#include <vector>
#include <cmath>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
using Eigen::Map;               	// 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
#include <RcppEigen.h>
using namespace Rcpp;
#define ZEPS 1e-10

// define simplex615
template <class F> // F is a function object
class simplex615 { // contains (dim+1) points of size (dim)

protected:
  std::vector<std::vector<double> > X; // (dim+1)*dim matrix
  std::vector<double> Y; // (dim+1) vector
  std::vector<double> midPoint; // variables for update
  std::vector<double> thruLine; // variables for update
  int dim, idxLo, idxHi, idxNextHi; // dimension, min, max, 2ndmax values
  void evaluateFunction(F& foo); // evaluate function value at each point
  void evaluateExtremes(); // determine the min, max, 2ndmax
  void prepareUpdate(); // calculate midPoint, thruLine
  bool updateSimplex(F& foo, double scale); // for reflection/expansion..
  void contractSimplex(F& foo); // for multiple contraction
  static int check_tol(double fmax, double fmin, double ftol); // check tolerance

  // attempting to give data to simplex
  Eigen::Map<Eigen::MatrixXd> &data;
  int column;
  std::vector<double> weights;


public:
  simplex615(double* p, int d, Eigen::Map<Eigen::MatrixXd> &data, int column, std::vector<double> weights); // constructor with initial points
  void amoeba(F& foo, double tol); // main function for optimization
  std::vector<double>& xmin(); // optimal x value
  double ymin(); // optimal y value
};

template <class F> simplex615<F>::simplex615(double* p, int d, Eigen::Map<Eigen::MatrixXd> &dat, int col, std::vector<double> wts) : dim(d), data(dat) { // set dimension
  // Determine the space required
  X.resize(dim+1); // X is vector-of-vector, like 2-D array
  Y.resize(dim+1); // Y is function value at each simplex point
  midPoint.resize(dim); thruLine.resize(dim);
  for(int i=0; i < dim+1; ++i) { X[i].resize(dim); // allocate the size of content in the 2-D array
  } // Initially, make every point in the simplex identical
  for(int i=0; i < dim+1; ++i)
    for(int j=0; j < dim; ++j) X[i][j] = p[j]; // set each simple point to the starting point // then increase each dimension by one unit except for the last point
  for(int i=0; i < dim; ++i) X[i][i] += 1.; // this will generate a simplex
  // attempting to give data to simplex
  // data = dat;
  column = col;
  weights = wts;
}

template <class F>
void simplex615<F>::evaluateFunction(F& foo) {
  for(int i=0; i < dim+1; ++i) {
    Y[i] = foo(X[i], data, column, weights); // foo is a function object, which will be visited later
  }
}

template <class F>
void simplex615<F>::evaluateExtremes() {
  if ( Y[0] > Y[1] ) { // compare the first two points
    idxHi = 0; idxLo = idxNextHi = 1;
  }
  else {
    idxHi = 1; idxLo = idxNextHi = 0;
  }
  // for each of the next points
  for(int i=2; i < dim+1; ++i) {
    if ( Y[i] <= Y[idxLo] ) // update the best point if lower
      idxLo = i;
    else if ( Y[i] > Y[idxHi] ) { // update the worst point if higher
      idxNextHi = idxHi; idxHi = i;
    }
    else if ( Y[i] > Y[idxNextHi] ) // update also if it is the 2nd-worst point
      idxNextHi = i;
  }
}

template <class F> void simplex615<F>::prepareUpdate() {
  for(int j=0; j < dim; ++j) {
    midPoint[j] = 0; // average of all points but the worst point
  }
  for(int i=0; i < dim+1; ++i) {
    if ( i != idxHi ) { // exclude the worst point
      for(int j=0; j < dim; ++j) {
        midPoint[j] += X[i][j]; }
    }
  }
  for(int j=0; j < dim; ++j) { midPoint[j] /= dim;
    thruLine[j] = X[idxHi][j] - midPoint[j]; // direction for optimization
  }
}

template <class F>
bool simplex615<F>::updateSimplex(F& foo, double scale) {
  std::vector<double> nextPoint; // next point to evaluate
  nextPoint.resize(dim);
  for(int i=0; i < dim; ++i) {
    nextPoint[i] = midPoint[i] + scale * thruLine[i];
  }
  double fNext = foo(nextPoint,data,column,weights);
  if ( fNext < Y[idxHi] ) { // update only maximum values (if possible)
    for(int i=0; i < dim; ++i) { // because the order can be changed with
      X[idxHi][i] = nextPoint[i]; // evaluateExtremes() later
    }
    Y[idxHi] = fNext;
    return true;
  }
  else {
    return false; // never mind if worse than the worst
  }
}

template <class F> void simplex615<F>::contractSimplex(F& foo) {
  for(int i=0; i < dim+1; ++i) { if ( i != idxLo ) { // except for the minimum point
    for(int j=0; j < dim; ++j) { X[i][j] = 0.5*( X[idxLo][j] + X[i][j] ); // move the point towards minimum
    }
    Y[i] = foo(X[i], data, column, weights); // re-evaluate the function
  }
  }
}


template <class F>
void simplex615<F>::amoeba(F& foo, double tol) {
  evaluateFunction(foo); // evaluate the function at the initial points
  while(true) {
    evaluateExtremes(); // determine three important points
    prepareUpdate(); // determine direction for optimization
    if ( check_tol(Y[idxHi],Y[idxLo],tol) ) break; // check convergence
    updateSimplex(foo, -1.0); // reflection
    if ( Y[idxHi] < Y[idxLo] ) {
      updateSimplex(foo, -2.0); // expansion
    }
    else if ( Y[idxHi] >= Y[idxNextHi] ) {
      if ( !updateSimplex(foo, 0.5) ) { // 1-d contraction
        contractSimplex(foo); // multiple contractions
      }
    }
  }
}
template <class F> int simplex615<F>::check_tol(double fmax, double fmin, double ftol) {
  // calculate the difference
  double delta = fabs(fmax - fmin); // calculate the relative tolerance
  double accuracy = (fabs(fmax) + fabs(fmin)) * ftol; // check if difference is within tolerance
  return (delta < (accuracy + ZEPS));
}

template <class F>
std::vector<double>& simplex615<F>::xmin(){
  return X[idxLo];
}

template <class F>
double simplex615<F>::ymin(){
  return Y[idxLo];
}

// likelihoodFunc
class likelihoodFunc {
public:
  double operator() (std::vector<double>& beta, Eigen::Map<Eigen::MatrixXd> &data, int column, std::vector<double> weights) {
    int y = data.cols()-1;
    double likelihood =  0;
    // data(i,y) is the same as saying "Y_i"
    // data(i, column) is the same as saying "X_i"; note the index starts at 0 and not 1 for column
    for (int i = 0; i < data.rows(); i++) {
      likelihood -= weights[i]*(data(i, y)*(beta[0]+beta[1]*data(i, column)-log(1+exp(beta[0]+beta[1]*data(i, column)))) - (1-data(i,y))*(log(1+exp(beta[0]+beta[1]*data(i, column)))));
    }
    return likelihood;
  }
};

// univLogReg
// [[Rcpp::export]]
std::vector<double> univLogReg(std::vector<double> weights, Eigen::Map<Eigen::MatrixXd> &data, int column) {
  // data.col(data.cols()-1) is the last column, which stores the dependent variable
  // data.col(column) is the univariate predictor

  double point[2] = {0, 0}; // initial point to start
  likelihoodFunc foo; // WILL BE DISCUSSED LATER
  simplex615<likelihoodFunc> simplex(point, 2, data, column, weights); // create a simplex
  simplex.amoeba(foo, 1e-7); // optimize for a function

  std::vector<double> beta;
  beta.push_back(simplex.xmin()[0]);
  beta.push_back(simplex.xmin()[1]);

  /*
  double prob;
  for (int i = 0; i < data.rows(); i++) {
    prob = 1/(1 + exp(-1*(beta[0] + beta[1]*data(i, column))));
    hypothesis.push_back(prob >= 0.5);
  } */

  return beta;
}


// adaboost
// [[Rcpp::export]]
NumericMatrix adaboost(Eigen::Map<Eigen::MatrixXd> &data) {

  // initialize weights
  std::vector<double> weights;
  double l = data.rows();
  double init_weight = 1/l;
  for (int i = 0; i < l; i++) {
    weights.push_back(init_weight);
  }

  // Do for t = 1,..., T
  double T = data.cols()-1;
  std::vector<double> all_bt;
  double sum_bt = 0;

  // initialize object to store hypothesis; consider using a matrix or anything less stupid than this
  std::vector< std::vector<bool> > all_hypothesis;
  // initialize object to store betas
  std::vector< std::vector<double> > all_betas;

  for (int t = 0; t < T; t++) {
    // train the classifier with respect to the weighted sample set and obtain hypothesis
    std::vector<double> beta = univLogReg(weights, data, t); // t is the column
    std::vector<bool> hypothesis;
    double prob;
    for (int i = 0; i < data.rows(); i++) {
      prob = 1/(1 + exp(-1*(beta[0] + beta[1]*data(i, t))));
      hypothesis.push_back(prob >= 0.5);
    }
    all_hypothesis.push_back(hypothesis);
    all_betas.push_back(beta);

    // calculate the training error of hypothesis
    double training_error = 0;
    bool incorrect = false;
    double delta = 0.5; // This is a given number that terminates loop if the training error is too large

    for (int i = 1; i < l; i++) {
      incorrect = (hypothesis[i] != data(i, T) );
      training_error += weights[i]*incorrect;
    }

    // Worst training error = 0.5"
     if (training_error >= delta) {
      for(int i = 0; i < hypothesis.size(); i++) {
        if (hypothesis[i] == false) {
          hypothesis[i] = true;
        } else {
          hypothesis[i] = false;
        }
      }
      training_error = 1 - training_error;
    }

    // set b_t
    double b_t = log( (1-training_error)/training_error );
    all_bt.push_back(b_t);
    sum_bt += fabs(b_t);

    // update weights
    bool correct = false;
    double sum_weights = 0;
    for (int i = 0; i < l; i++) {
      correct = (hypothesis[i] == data(i,T) );
      weights[i] = weights[i]*exp(-1*b_t*correct);
      sum_weights += weights[i];
    }
    // scale the weights; we are looping through again which might be slow? consider using "transform" function instead
    for (int i = 0; i < l; i++) {
      weights[i] = weights[i]/sum_weights;
    }
  }

  // create vector c_ts and calculate final labels; store c_ts in the all_bt vector
  std::vector<double> labels;

  for (int i = 0; i < all_bt.size(); i++) {
    all_bt[i] = all_bt[i]/sum_bt;
  }

  // calculate final probabilities and assign labels
  double sum = 0;
  std::vector<bool> f_x; // this stores the final labels
  for (int i = 0; i <= l; i++) {
    sum = 0;
    for (int t = 0; t < T; t++) {
      sum += all_bt[t]*all_hypothesis[t][i];
    }
    f_x.push_back(sum >= 0.5);
  }
  // Note: what exactly should we be outputting?
  // Once we hook this up to R, we should be outputting a trained classifier that can then
  // be used for predicted; e.g. we should be outputting the betas from the logistic regression
  // and all_bt (which is c_t)
  // For each new point we want to predict, we calculate the hypothesis at that point using all sets of
  // beta values (e.g. if we have 6 X columns we end up with 6 hypotheses at that point). Then, our c_t
  // tells us which hypthosis to use at that point.
  // Note that this means that our univLogReg function should probably be returning the beta's, not the
  // explicit hypotheses.

  // initialize object to store betas and c_t's
  NumericMatrix results(T,3);
  // std::cout << "Please God" << std::endl << "T: " << T << std::endl;
  for (int i = 0; i < T; i++) {
    results(i,0) = all_betas[i][0];
    results(i,1) = all_betas[i][1];
    results(i,2) = all_bt[i];
  }
  return results;
}


// predictAda615
// [[Rcpp::export]]
NumericVector predictAda615(NumericMatrix test, NumericMatrix train_results) {
  NumericMatrix hypothesis(test.nrow(), train_results.nrow());
  double prob;
  NumericVector predicted(test.nrow());
  for (int i = 0; i < test.nrow(); i++) {
    for (int j = 0; j < train_results.nrow(); j++) {
      prob = 1/(1 + exp(-1*(train_results(j,0) + train_results(j,1)*test(i,j))));
      // std::cout << "prob: " << j << ": "<< prob << std::endl;
      if (prob >= 0.5) {
        predicted(i) += 1*train_results(j,2);
      }
      if (predicted(i) >= 0.5) {
        predicted(i) = 1;
      } else {
        predicted(i) = 0;
      }
    }
  }
  return predicted;
}


// error rate
// [[Rcpp::export]]
double errorRate(NumericVector predicted, NumericVector label){
  if (predicted.size() != label.size()) {
    std::cout << "Error: size of predicted and label vectors do not match.";
    return 0;
  }
  double error=0;
  for(int i = 0;i < label.size(); i++) {
    if (predicted[i] != label[i]) {
      error += 1;
    }
  }
  return error/label.size();
}
