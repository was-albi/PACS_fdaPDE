#ifndef FPIRLS_IMP_H
#define FPIRLS_IMP_H

#include "FPIRLS.h"

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim>& mesh, InputHandler& inputData, VectorXr mu0):mesh_(mesh), inputData_(inputData){
  Rprintf("Hello I'm FPIRLS_Base constructor \n");

  //mu_(inputData.getLambda().size(),mu0);
  //current_J_values(inputData.getLambda().size(),std::array<Real,2>{1,1});
  //past_J_values(inputData.getLambda().size(),std::array<Real,2>{1,1});

  std::string saving_filename = "TEST_COSTRUTTORE";
  saving_filename = saving_filename + ".txt";
  printer::saveVectorXr(saving_filename,mu0);


  //MU,J
  for(UInt j=0; j< inputData.getLambda().size() ; j++){
    mu_.push_back(mu0);
    current_J_values.push_back(std::array<Real,2>{1,1});
    past_J_values.push_back(std::array<Real,2>{1,1});
  }


}; // Constructor


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::apply( const ForcingTerm& u){
  // f-PRILS implementation

  this->inputData_.copyInitialObservations();
  this->inputData_.copyLambdaVector();
  this->inputData_.setDOFflag();

  // STEP 0: initialization
  //  if(mu_.size()==0) //mu can be initialized or not by the user
      //initialize_mu(inputData_.getObservations());

  G_.resize(mu_.size());
  WeightsMatrix_.resize(mu_.size());
  _beta_hat.resize(mu_.size());
  _fn_hat.resize(mu_.size());
  _dof.resize(mu_.size());
  _solution.resize(mu_.size());
  pseudoObservations_.resize(mu_.size());
  n_iterations = std::vector<UInt>(mu_.size(),0);
  _GCV.resize(mu_.size(),-1);

  FiniteElement<Integrator, ORDER, mydim, ndim> fe;

  typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
  Assembler::operKernel(mass, mesh_, fe, R0_);

  if(isSpaceVarying)
  {
  	Assembler::forcingTerm(mesh_, fe, u, forcingTerm);
  }

  std::string saving_filename = "TEST_1";
  saving_filename = saving_filename + ".txt";
  printer::saveVectorXr(saving_filename,mu_[0]);

for(UInt i=0 ; i < mu_.size() ; i++){

  current_J_values[i][0] = past_J_values[i][0] + 2*inputData_.get_treshold();
  current_J_values[i][1] = past_J_values[i][1] + 2*inputData_.get_treshold();


    saving_filename = "TEST_2";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);



  this->inputData_.setCurrentLambda(i); // set right lambda for the iteration

  while(stopping_criterion(i))  // n_iterations < inputData_.get_maxiter() && past_J_value - current_J_value < inputData_.get_treshold()
  {

  Rprintf("FPIRLS.apply: while iteration number: %d \n", n_iterations[i]);
  // STEP (1)

  saving_filename = "TEST_3";
  saving_filename = saving_filename + ".txt";
  printer::saveVectorXr(saving_filename,mu_[0]);


    compute_G(i);

    saving_filename = "TEST_4";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);


    compute_Weights(i);

    saving_filename = "TEST_5";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);


    compute_pseudoObs(i);


    saving_filename = "TEST_6";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);


  // STEP (2)
    this->inputData_.updatePseudodata(pseudoObservations_[i], WeightsMatrix_[i]);
    update_solution(i); // here I'm performing FERegression using appropriate objects. I need pseudo data and mesh and it computes solution and dof


    saving_filename = "TEST_7";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);

  // STEP (3)

    compute_mu(i);


    saving_filename = "TEST_8";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);

  // update J

    past_J_values[i] = current_J_values[i];
    current_J_values[i] = compute_J(i);

    n_iterations[i]++;

    saving_filename = "TEST_9";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);


  } //end while

  _J_minima.push_back(current_J_values[i][0]+current_J_values[i][1]);

  if(this->inputData_.getDOF_GAM()){
    saving_filename = "TEST_GCV";
    saving_filename = saving_filename + ".txt";
    printer::saveVectorXr(saving_filename,mu_[0]);
    compute_GCV(i);
  }

}// end for


saving_filename = "TEST_10";
saving_filename = saving_filename + ".txt";
printer::saveVectorXr(saving_filename,mu_[0]);


}// end apply


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::update_solution(UInt& lambda_index){
  // performs step (2) of PIRLS. It requires pseudo data after step(1) and mimic regression skeleton behaviour

  // Here we have to solve a weighted regression problem with laplacian penality
  // COMPOSE PSEUDO DATA
  // I SetMethods di Regression data sono protected, perciò non posso modificare l'oggetto. Tuttavia, potrei creare un metodo ad-hoc in RegressiondDataGAM che faccia ciò di cui ho bisogno, ovvero:
  //      1. l'oggetto di tipo Regression Data abbia tutti i parametri necessari per compiere la regressione
  //      2. observations substituted by z, WeightsMatrix substituted by W

  //Set up regression
  MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim> regression(mesh_, inputData_);

  regression.apply();

  _solution[lambda_index] = regression.getSolution()[0];
  _dof[lambda_index] = regression.getDOF()[0];

  std::string saving_filename = "solution_entire_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,_solution[lambda_index]);

/*  VectorXr solution_coefs = VectorXr::Zero(_solution[lambda_index].size()/2);
  VectorXr laplacian_coefs = VectorXr::Zero(_solution[lambda_index].size()/2);

  for(UInt i=0; i < _solution[lambda_index].size()/2; i++){
    solution_coefs[i] = _solution[lambda_index][i];
    laplacian_coefs[i] = _solution[lambda_index][i + _solution[0].size()/2];
  }


  		saving_filename = "apply_return_solution_coefs_";
  		saving_filename = saving_filename + ".txt";
  		printer::saveVectorXr(saving_filename,solution_coefs);

  		saving_filename = "apply_return_laplacian_coefs_";
  	  saving_filename = saving_filename + ".txt";
  	  printer::saveVectorXr(saving_filename,laplacian_coefs);*/

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_pseudoObs(UInt& lambda_index){

  VectorXr g_mu;
  VectorXr first_addendum;
  VectorXr y = inputData_.getInitialObservations();

  g_mu.resize(mu_[lambda_index].size());
  first_addendum.resize(mu_[lambda_index].size());

  for(auto i=0; i < mu_[lambda_index].size(); i++){
    g_mu(i) = link(mu_[lambda_index](i));
    first_addendum(i) = G_[lambda_index](i)*(y(i)-mu_[lambda_index](i));
  }

  pseudoObservations_[lambda_index] = first_addendum + g_mu;

  std::string saving_filename = "pseudoObservations_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,pseudoObservations_[lambda_index]);
}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_G(UInt& lambda_index){

  G_[lambda_index].resize(mu_[lambda_index].size());

  for(UInt i = 0; i<mu_[lambda_index].size(); i++){
    G_[lambda_index](i) = link_deriv(mu_[lambda_index](i));
  }

  std::string saving_filename = "G_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,G_[lambda_index]);

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_Weights(UInt& lambda_index){
// computed W elementwise (it is a diagonal matrix)

  WeightsMatrix_[lambda_index].resize( mu_[lambda_index].size());

  for(auto i=0; i < mu_[lambda_index].size(); i++){
    WeightsMatrix_[lambda_index](i) = 1/(pow(G_[lambda_index](i),2)*(var_function( mu_[lambda_index](i))));
  }

  std::string saving_filename = "WeightsMatrix_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,WeightsMatrix_[lambda_index]);

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_mu(UInt& lambda_index){


  VectorXr X_beta = VectorXr::Zero(mu_[lambda_index].size());

  // METODO STATIC: OLD VERSION
  //MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>::computeRegressionSolution(_solution, pseudoData, mesh_, fn_hat, beta_hat);

  RegressionEstimates<InputHandler, ORDER, mydim, ndim> regressionEst(inputData_, mesh_, _solution[lambda_index]);

  regressionEst.computeEstimates();

  _beta_hat[lambda_index] = regressionEst.getBetaEst()[0];
  _fn_hat[lambda_index] = regressionEst.getFunctionEst()[0];

  std::string saving_filename;

  saving_filename = "compute_mu_betahat_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,_beta_hat[lambda_index]);

  saving_filename = "compute_mu_fnhat_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,_fn_hat[lambda_index]);

  /*  for(UInt i = 0, i < pseudoData.getLambda().size(), i++){

        X_beta = inputData_.getCovariates()*beta_hat[i];

        for(UInt j=0; j < X_beta.size() ; j++){

          mu_(j) = inv_link(X_beta[j] + fn_hat[i][j]);

        }// end for-j

      }// end for-i


  */

  if(inputData_.getCovariates().rows()>0)
    X_beta = inputData_.getCovariates()*_beta_hat[lambda_index];

  saving_filename = "X_beta";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,X_beta);


  for(UInt j=0; j < X_beta.size() ; j++){
      mu_[lambda_index](j) = inv_link(X_beta[j] + _fn_hat[lambda_index][j]);
  }

  saving_filename = "mu_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,mu_[lambda_index]);

}//end method


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
bool FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::stopping_criterion(UInt& lambda_index){
  // return true if the f-PIRLS has to perform another iteration, false if it has to be stopped

  bool do_stop_by_iteration = false;  // Do I need to stop becouse n_it > n_max?
  bool do_stop_by_treshold = false; // Do I need to stop becouse |J(k) - J(k+1)| < treshold?

  if(n_iterations[lambda_index] > inputData_.get_maxiter()){
    do_stop_by_iteration = true;
  }

  if(n_iterations[lambda_index] > 1){
    if(abs(past_J_values[lambda_index][0]+past_J_values[lambda_index][1] - current_J_values[lambda_index][0] - current_J_values[lambda_index][1]) < inputData_.get_treshold()){
        do_stop_by_treshold = true;
   }
  }


  return !(do_stop_by_iteration || do_stop_by_treshold );
}



// Laplace or Elliptic
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<InputHandler,Integrator,ORDER, mydim, ndim>::apply(){

  FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::apply(ForcingTerm(std::vector<Real>(1)));

}

// SpaceVarying
template <typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS<GAMDataEllipticSpaceVarying,Integrator,ORDER, mydim, ndim>::apply(){

  this->isSpaceVarying = true;
  FPIRLS_Base<GAMDataEllipticSpaceVarying,Integrator,ORDER, mydim, ndim>::apply(this->inputData_.getU());

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Scaled<InputHandler, Integrator, ORDER, mydim, ndim>::apply(){

  FPIRLS<InputHandler, Integrator, ORDER, mydim, ndim>::apply();

  // Phi estimate
  if(scale_parameter_flag_ && this->_dof[0]>=0 ){ //if scale_flag is true and the dofs have been estimated
    compute_scale_param();
  }

}

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Scaled<InputHandler, Integrator, ORDER, mydim, ndim>::compute_scale_param(){

  _scale_parameter_estimates.resize(this->_dof.size());
  const UInt n_obs = this->inputData_.getObservations().size();

  for(UInt i=0; i < this->_dof.size();i++){
    _scale_parameter_estimates[i] = this->current_J_values[i][0]/(n_obs - this->_dof[i]);
  }

}

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
std::array<Real,2> FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_J(UInt& lambda_index){

  Real parametric_value = 0;
  Real non_parametric_value = 0;
  Real tmp;


  VectorXr Lf;

  const VectorXr y = inputData_.getInitialObservations();

  for(UInt i=0; i < mu_.size(); i++){
    tmp = sqrt( var_function( mu_[lambda_index](i)) ) * (y(i) - mu_[lambda_index](i)) ;
    parametric_value += tmp*tmp;
  }
  Rprintf("\t \t norm_value: %f \n", parametric_value);
  // not sure if the following work/is correct


  std::string saving_filename = "R0_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::SaveMatrixXr(saving_filename,R0_);

  Lf.resize(_solution[lambda_index].size()/2);
  for(UInt i=0; i< Lf.size(); i++){
    Lf(i) = _solution[lambda_index](Lf.size() + i);
  }

  saving_filename = "Lf_";
  saving_filename = saving_filename + std::to_string(lambda_index) + ".txt";
  printer::saveVectorXr(saving_filename,Lf);


  if(isSpaceVarying)
  {
      Lf = Lf - forcingTerm;
  }

  non_parametric_value = Lf.transpose() * R0_ * Lf;
  non_parametric_value = inputData_.getGlobalLambda()[lambda_index]*non_parametric_value;

  std::array<Real,2> returnObject{parametric_value, non_parametric_value};

  return returnObject;
}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_GCV(UInt& lambda_index){

//GCV COMPUTATION
MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim> regression(mesh_, inputData_);

regression.computeDegreesOfFreedom();
_dof[lambda_index] = regression.getDOF()[0];

std::string saving_filename = "DIMENSIONE_dof_computeGCV";
saving_filename = saving_filename + ".txt";
printer::SaveDimension(saving_filename,_dof);

VectorXr y = inputData_.getInitialObservations();
Real GCV_value = 0;

for(UInt j=0; j < y.size();j++)
  GCV_value += dev_function(mu_[lambda_index][j], y[j]); //norm computation

  saving_filename = "DIMENSIONE_GCV_value";
  saving_filename = saving_filename + ".txt";
  printer::SaveDimension(saving_filename,_dof);

GCV_value *= y.size();

GCV_value /= (y.size()-inputData_.getTuneParam()*_dof[lambda_index])*(y.size()-inputData_.getTuneParam()*_dof[lambda_index]);

saving_filename = "DIMENSIONE_GCV_value_2";
saving_filename = saving_filename + ".txt";
printer::SaveDimension(saving_filename,_dof);

_GCV[lambda_index] = GCV_value;

saving_filename = "DIMENSIONE_GCV_value_3";
saving_filename = saving_filename + ".txt";
printer::SaveDimension(saving_filename,_dof);

}

/*
template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Bernoulli<InputHandler,Integrator,ORDER, mydim, ndim>::initialize_mu(const VectorXr& y){
      this->mu_.resize(y.size());
      for(UInt i = 0; i < y.size(); i++){
          this->mu_[i] = 0.5 * (y[i] + 0.5); // It is different for binary or non-binary outcomes
      }
}
*/




#endif
