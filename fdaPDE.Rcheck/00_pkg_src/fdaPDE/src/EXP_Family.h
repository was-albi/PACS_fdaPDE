#ifndef EXP_FAMILY_H
#define EXP_FAMILY_H

#include <cmath>
#include <math.h>
#include "GLM.h"

template <typename Data_Type, typename Parameter_Type>
class Exp_Family {

protected:

  double scale_parameter = 1; // default value since it is 1 for most of the distributions
  int scale_parameter_flag = 0; // It tells me if I have to estimate phi (W. do it only for gamma)

public:

  virtual Parameter_Type link(Parameter_Type& mu) = 0; // g(.)

  virtual Parameter_Type link_deriv(Parameter_Type& mu) = 0; // g'(.)

  virtual Parameter_Type inv_link(Parameter_Type& theta) = 0; // g^-1(.)

  virtual Parameter_Type var_function(Parameter_Type& mu) = 0; // V(mu)

};

// Bernoulli

class Bernoulli_Dist : public Exp_Family<Eigen::VectorXi, double > {

public:

  double link(double& mu){

      double theta = log(mu/(1 - mu));
      return theta;

  };

  double inv_link(double& theta){

    double mu = exp(theta)/(1 + exp(theta));
    return mu;

  };

  double link_deriv(double& mu){ return 1/(mu*(1-mu)); };

  double var_function(double& mu){ return(mu*(1-mu)); };

};


// Poisson

class Poisson_Dist : public Exp_Family<Eigen::VectorXi, double> {

public:

  double link(double& mu){
      double theta = log(mu);
      return theta;
  };

  double link_deriv(double& mu){ return 1/mu; };

  double inv_link(double& theta){
    double mu = exp(theta);
    return mu;
  }


  double var_function(double& mu){ return mu ;};

};

// Exponential

class Exp_Dist : public Exp_Family<Eigen::VectorXd, double> {

public:

  double link(double& mu){ return -1/mu ; };

  double link_deriv(double& mu){ return 1/(mu*mu); };

  double inv_link(double& theta){ return -1/mu; }

  double var_function(double& mu){ return mu*mu ;};

};

// Gamma

class Gamma_Dist : public Exp_Family<Eigen::VectorXd, double> {

public:

  Gamma_Dist(double scale_parameter_flag_ , double scale_parameter_){     // w. treats differently gamma scale parameter from others
    scale_parameter = scale_parameter_;
    scale_parameter_flag = scale_parameter_flag_;
  };

  double link(double& mu){ return -1/mu ; };

  double link_deriv(double& mu){ return 1/(mu*mu); };

  double inv_link(double& theta){ return -1/mu; }

  double var_function(double& mu){ return mu*mu ;};

};


// Probit
// TO BE DONE: how to get normal cdf and its inverse? thay're needed for probit

// Cloglog
// TO BE DONE: W. write his wrong in the old R code, has to be fixed here

#endif
