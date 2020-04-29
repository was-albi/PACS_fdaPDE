#ifndef FPIRLS_IMP_H
#define FPIRLS_IMP_H

#include "FPIRLS.h"


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::apply( const ForcingTerm& u){
  // f-PRILS implementation

  this->inputData_.copyInitialObservations();

  // STEP 0: initialization
  Rprintf("FPIRLS.apply: STEP 0: initialization \n");
  //  if(mu_.size()==0) //mu can be initialized or not by the user
      //initialize_mu(inputData_.getObservations());

  G_.resize(mu_.size());
  WeightsMatrix_.resize(mu_.size());
  n_iterations = 0;
  current_J_value[0] = past_J_value[0] + 2*inputData_.get_treshold();
  current_J_value[1] = past_J_value[1] + 2*inputData_.get_treshold();


  FiniteElement<Integrator, ORDER, mydim, ndim> fe;

  typedef EOExpr<Mass> ETMass; Mass EMass; ETMass mass(EMass);
  Assembler::operKernel(mass, mesh_, fe, R0_);

  if(isSpaceVarying)
  {
  	Assembler::forcingTerm(mesh_, fe, u, forcingTerm);
  }

  while(stopping_criterion())  // n_iterations < inputData_.get_maxiter() && past_J_value - current_J_value < inputData_.get_treshold()
  {

  Rprintf("FPIRLS.apply: while iteration number: %d \n", n_iterations);
  // STEP (1)

    compute_G();

    compute_Weights();

    compute_pseudoObs();

    Rprintf("\t \t mu_(0): %f, size: %d \n", mu_(0), mu_.size() );

  // STEP (2)
    this->inputData_.updatePseudodata(pseudoObservations_, WeightsMatrix_);
    update_solution(); // here I'm performing FERegression using appropriate objects. I need pseudo data and mesh and it computes solution and dof

  // STEP (3)

    compute_mu();

  // update J

    past_J_value = current_J_value;
    current_J_value = compute_J();
    Rprintf("\t \t current_J_value: param -> %f non-param -> %f \n", current_J_value[0],current_J_value[1]);

    n_iterations++;

  } //end while


  _J_minima.push_back(current_J_value[0]+current_J_value[1]);

}// end apply


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::update_solution(){
  // performs step (2) of PIRLS. It requires pseudo data after step(1) and mimic regression skeleton behaviour

  // Here we have to solve a weighted regression problem with laplacian penality
  // COMPOSE PSEUDO DATA
  // I SetMethods di Regression data sono protected, perciò non posso modificare l'oggetto. Tuttavia, potrei creare un metodo ad-hoc in RegressiondDataGAM che faccia ciò di cui ho bisogno, ovvero:
  //      1. l'oggetto di tipo Regression Data abbia tutti i parametri necessari per compiere la regressione
  //      2. observations substituted by z, WeightsMatrix substituted by W

  //Set up regression
  MixedFERegression<InputHandler, Integrator, ORDER, mydim, ndim> regression(mesh_, inputData_);

  regression.apply();
  //  regression.apply();
  _solution = regression.getSolution();
  _dof = regression.getDOF();

  std::string saving_filename = "solution_entire_";
  saving_filename = saving_filename + std::to_string(n_iterations) + ".txt";
  printer::saveVectorXr(saving_filename,_solution[0]);

  VectorXr solution_coefs = VectorXr::Zero(_solution[0].size()/2);
  VectorXr laplacian_coefs = VectorXr::Zero(_solution[0].size()/2);

  for(UInt i=0; i < _solution[0].size()/2; i++){
    solution_coefs[i] = _solution[0][i];
    laplacian_coefs[i] = _solution[0][i + _solution[0].size()/2];
  }


  		saving_filename = "apply_return_solution_coefs_";
  		saving_filename = saving_filename + ".txt";
  		printer::saveVectorXr(saving_filename,solution_coefs);

  		saving_filename = "apply_return_laplacian_coefs_";
  	  saving_filename = saving_filename + ".txt";
  	  printer::saveVectorXr(saving_filename,laplacian_coefs);

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_pseudoObs(){

  VectorXr g_mu;
  VectorXr first_addendum;
  VectorXr y = inputData_.getInitialObservations();

  g_mu.resize(mu_.size());
  first_addendum.resize(mu_.size());

  for(auto i=0; i < mu_.size(); i++){
    g_mu(i) = link(mu_(i));
    first_addendum(i) = G_(i)*(y(i)-mu_(i));
  }

  pseudoObservations_ = first_addendum + g_mu;

  std::string saving_filename = "pseudoObservations_";
  saving_filename = saving_filename + std::to_string(n_iterations) + ".txt";
  printer::saveVectorXr(saving_filename,pseudoObservations_);
}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_G(){

  for(UInt i = 0; i<mu_.size(); i++){
    G_(i) = link_deriv(mu_(i));
  }

  std::string saving_filename = "G_";
  saving_filename = saving_filename + std::to_string(n_iterations) + ".txt";
  printer::saveVectorXr(saving_filename,G_);

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_Weights(){
// computed W elementwise (it is a diagonal matrix)

  for(auto i=0; i < mu_.size(); i++){
    WeightsMatrix_(i) = 1/(pow(G_(i),2)*(var_function( mu_(i))));
  }

  std::string saving_filename = "WeightsMatrix_";
  saving_filename = saving_filename + std::to_string(n_iterations) + ".txt";
  printer::saveVectorXr(saving_filename,WeightsMatrix_);

}


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
void FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_mu(){


  VectorXr X_beta = VectorXr::Zero(mu_.size());

  // METODO STATIC: OLD VERSION
  //MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>::computeRegressionSolution(_solution, pseudoData, mesh_, fn_hat, beta_hat);

  RegressionEstimates<InputHandler, ORDER, mydim, ndim> regressionEst(inputData_, mesh_, _solution);

  regressionEst.computeEstimates();

  _beta_hat = regressionEst.getBetaEst();
  _fn_hat = regressionEst.getFunctionEst();

  std::string saving_filename;

  saving_filename = "compute_mu_betahat_";
  saving_filename = saving_filename + std::to_string(n_iterations) + ".txt";
  printer::saveVectorXr(saving_filename,_beta_hat[0]);

  saving_filename = "compute_mu_fnhat_";
  saving_filename = saving_filename + std::to_string(n_iterations) + ".txt";
  printer::saveVectorXr(saving_filename,_fn_hat[0]);

  /*  for(UInt i = 0, i < pseudoData.getLambda().size(), i++){

        X_beta = inputData_.getCovariates()*beta_hat[i];

        for(UInt j=0; j < X_beta.size() ; j++){

          mu_(j) = inv_link(X_beta[j] + fn_hat[i][j]);

        }// end for-j

      }// end for-i


  */

  if(inputData_.getCovariates().rows()>0)
    X_beta = inputData_.getCovariates()*_beta_hat[0];

  saving_filename = "X_beta";
  saving_filename = saving_filename + std::to_string(n_iterations) + ".txt";
  printer::saveVectorXr(saving_filename,X_beta);


  for(UInt j=0; j < X_beta.size() ; j++){
      mu_(j) = inv_link(X_beta[j] + _fn_hat[0][j]);
  }

  saving_filename = "mu_";
  saving_filename = saving_filename + std::to_string(n_iterations) + ".txt";
  printer::saveVectorXr(saving_filename,mu_);

}//end method


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
bool FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::stopping_criterion(){
  // return true if the f-PIRLS has to perform another iteration, false if it has to be stopped

  bool do_stop_by_iteration = false;  // Do I need to stop becouse n_it > n_max?
  bool do_stop_by_treshold = false; // Do I need to stop becouse |J(k) - J(k+1)| < treshold?

  if(n_iterations > inputData_.get_maxiter()){
    do_stop_by_iteration = true;
  }

  if(n_iterations > 1){
    if(abs(past_J_value[0]+past_J_value[1] - current_J_value[0] - current_J_value[1]) < inputData_.get_treshold()){
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
    _scale_parameter_estimates[i] = this->current_J_value[0]/(n_obs - this->_dof[i]);
  }

}

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
std::array<Real,2> FPIRLS_Base<InputHandler,Integrator,ORDER, mydim, ndim>::compute_J(){

  Real parametric_value = 0;
  Real non_parametric_value = 0;
  Real tmp;


  VectorXr Lf;

  const VectorXr y = inputData_.getInitialObservations();

  for(UInt i=0; i < mu_.size(); i++){
    tmp = sqrt( var_function( mu_(i)) ) * (y(i) - mu_(i)) ;
    parametric_value += tmp*tmp;
  }
  Rprintf("\t \t norm_value: %f \n", parametric_value);
  // not sure if the following work/is correct


  std::string saving_filename = "R0_";
  saving_filename = saving_filename + std::to_string(n_iterations) + ".txt";
  printer::SaveMatrixXr(saving_filename,R0_);

  Lf.resize(_solution[0].size()/2);
  for(UInt i=0; i< Lf.size(); i++){
    Lf(i) = _solution[0](Lf.size() + i);
  }

  saving_filename = "Lf_";
  saving_filename = saving_filename + std::to_string(n_iterations) + ".txt";
  printer::saveVectorXr(saving_filename,Lf);


  if(isSpaceVarying)
  {
      Lf = Lf - forcingTerm;
  }

  non_parametric_value = Lf.transpose() * R0_ * Lf;
  non_parametric_value = inputData_.getLambda()[0]*non_parametric_value;

  std::array<Real,2> returnObject{parametric_value, non_parametric_value};

  return returnObject;
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
