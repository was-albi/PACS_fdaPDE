#ifndef FPIRLS_IMP_H
#define FPIRLS_IMP_H

#include <cmath>
#include <math.h>
#include "GLM.h"
#include "mixedFERegression.h"
#include "fdaPDE.h"
#include "FPIRLS.h"



template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::apply(){
// f-PRILS implementation

// STEP 0: initialization
initialize_mu(inputData_.getObservations());
G_.resize(mu_.size());
W_.resize(mu_.size());
pseudo_Obs.resize(mu_.size());
pseudo_Cov.resize(inputData_.getCovariates().size());
n_iterations_ = 0;
past_J_value = 1;
current_J_value = past_J_value + 2*inputData_.get_treshold();

while(stopping_criterion())  // n_iterations < inputData_.get_maxiter() && past_J_value - current_J_value < inputData_.get_treshold()
{

// STEP (1)
  compute_G();
  compute_W();
  compute_z();

// STEP (2)
// Here we have to solve a weighted regression problem with laplacian penality
// COMPOSE PSEUDO DATA
// I SetMethods di Regression data sono protected, perciò non posso modificare l'oggetto. Tuttavia, potrei creare un metodo ad-hoc in RegressiondDataGAM che faccia ciò di cui ho bisogno, ovvero:
//      1. l'oggetto di tipo Regression Data abbia tutti i parametri necessari per compiere la regressione
//      2. observation substituted with pseudo obs and multiplied by W^1/2
//      3. observation and design matrix multiplied by W^1/2


// !!! INIZIALIZE (CREATE ATM) PSEUDODATA HERE
for(auto i=0;i<inputData_.getObservations().size(); i++){
  pseudo_Obs[i] = root_W_[i]*inputData_.getObservations()[i];
  for(auto j=0; j<inputData_.getCovariates().cols(); j++){
    pseudo_Cov[i,j] = root_W_[i]*inputData_.getCovariates()[i,j];
  }
}

RegressionData pseudoData_(inputData_.getLocations(), pseudo_Obs,  inputData_.getOrder(),  inputData_.getLambda(), pseudo_Cov, inputData_.getIncidenceMatrix(), inputData_.getDirichletIndices(), inputData_.getDirichletValues() , inputData_.computeDOF());

update_solution(); // here I'm performing FERegression using appropriate objects. I need pseudo data and mesh and it computes solution and dof

// STEP (3)
compute_mu();

// update J
past_J_value = current_J_value;
current_J_value = compute_J();
n_iterations++;

} //end while

_J_minima.push_back(current_J_value);

}// end apply

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::update_solution(){
// performs step (2) of PIRLS. It requires pseudo data after step(1) and mimic regression skeleton behaviour

MixedFERegression<InputHandler, Integrator,ORDER, mydim, ndim> regression(mesh,pseudoData_);

regression.apply();
_solution = regression.getSolution();
_dof = regression.getDOF();

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_z(){

vectorXr g_mu;
vectorXr first_addendum;
vectorXr y = inputData_.getObservations();

g_mu.resize(mu_.size());
first_addendum.resize(mu_.size());

for(auto i=0; i < mu_.size(); i++){
  g_mu[i] = link(mu_[i]);
  first_addendum[i] = G_[i]*(y[i]-mu_[i]);
}

z_ = first_addendum + g_mu;

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_G(){

for(auto i = 0; i<mu_.size(); i++){
  G_[i] = link_deriv(mu_[i]);
}

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_W(){
// computed W^1/2 directly

for(auto i=0; i < mu_.size(); i++){
  root_W_[i] = 1/(G_[i]*math::sqrt(var_function(mu_[i])));
}

}



template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_mu(){

VextorXr f_evaluated, beta_coef, X_beta;
UInt Length = inputData_.getLocations().size();
std::vector<Real> X,Y; //devo settare la dimensione?
std::vector<boolean> isinside[Length];
VectorXr solution_coeff;

f_evaluated.resize(Length);
solution_coeff.resize(_solution[0].size()/2);

for(auto i = 0; i < Length; i++){
  X.push_back(inputData_.getLocations()[i][1]);
  Y.push_back(inputData_.getLocations()[i][2]);
}

for(auto i = 0 , i < _solution[0].size()/2 , i++){
 soution_coeff[i] = _solution[0][i];
} // take only the value of f and not of laplacian(f)

Evaluator<ORDER> f_evaluator(mesh_);

f_evaluator.eval(X[0],Y[0], Length, solution_coeff, true, f_evaluated[0], isinside[0]);

// solve weighted least squares
beta_coef = (root_W_*inputData_.getCovariates()).bdcSvd(ComputeThinU | ComputeThinV).solve(root_W_*(z_ - f_evaluated)) ;

X_beta = inputData_.getCovariates()*beta_coef;

for(auto i=0; i < Length ; i++){

  mu_[i] = inv_link(X_beta[i] + f_evaluated[i]);

}

}



template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
bool FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::stopping_criterion(){
// return true if the f-PIRLS has to perform another iteration, false if it has to be stopped

bool do_stop_by_iteration = false;  // Do I need to stop becouse n_it > n_max?
bool do_stop_by_treshold = false; // Do I need to stop becouse |J(k) - J(k+1)| < treshold?

if(n_iterations > inputData_.get_maxiter()){
  do_stop_by_iteration = true;
}

past_J_value = current_J_value;
current_J_value = compute_J();

if(n_iterations > 2)
 if(abs(past_J_value - current_J_value) < inputData_.get_treshold()){
  do_stop_by_treshold = true;
  }
}

return !(do_stop_by_iteration || do_stop_by_treshold);

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
Real FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_J(){

  Real norm_value = 0;
  Real tmp;

  Real m_value;
  FiniteElement<Integrator, ORDER, mydim, ndim> fe;
  SpMat R0_;
  VectorXr lapl_f;

  for(auto i=0; i < mu_.size() ; i++){

    tmp = math::sqrt(var_function(mu_[i])) * (inputData_.getObservations()[i] - mu_[i]);
    norm_value += tmp*tmp;
  }

  // not sure if the following work/is correct

  typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
  Assembler::operKernel(mass, mesh_, fe, R0_);

  lapl_f.resize(_solution[0].size()/2);
  for(auto i=0, i< lapl_f.size(), i++){
    lapl_f[i] = _solution[0][lapl_f.size() + i];
  }

  m_value = lapl_f.dot((R0_*lapl_f)) ;

  return (norm_value + inputData_.getLambda()[0]*m_value);

}
