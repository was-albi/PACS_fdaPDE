#ifndef _FPIRLS_H
#define _FPIRLS_H

#include <cmath>
#include <math.h>
#include "mixedFERegression.h"
#include "evaluator.h"
#include "fdaPDE.h"


template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_Base {

  protected:

  const MeshHandler<ORDER, mydim, ndim> &mesh_;
  const InputHandler& inputData_; // it contains the original data of the problem

  //RegressionData& pseudoData_; // contains data used to perform step (2) of f-PIRLS algo
  // question is: should I introduce a new template for this (data in step (2))? At the moment it is RegressionData, but one can imagine future developements around childs of regressionData...

  VectorXr mu_; // mean vector
  VectorXr z_; // pseudo observations
  VectorXr G_; // diag(link_deriv(mu)) it is a vector since it would be more memory consuming to declere it as a matrix
  VectorXr W_; // (G^-2 * Var_function(mu)^-1) it is a vector for the same reason as above

  Real current_J_value;
  Real past_J_value; // stores the value of the functional J at each iteration in order to apply the stopping criterion
  UInt n_iterations; // current n° of iteration of PIRLS
  //UInt max_num_iterations; // max n° of iterations allowed
  //double threshold; // limit in difference among J_k and J_k+1 for which we stop PIRLS algo
  // NOTA IMPORTANTE: le ultime 2 variabili andrebbero inserite in una classe figlia di regression data. Se così fosse vorrebbe dire che pseudoData e inputData non sono due oggetti della stessa classe. Infatti pseudoData dovrebbe potersi inserire all'interno di mixedFERegression"qualcosa" per effettuare lo step (2) dell'algo. In tal caso sarebbe molto comodo un casting tra la classe di inputData e quella di pseudoData che ignori i parametri aggiuntivi e permetta di riciclare gli algo di mixedFERegression.
  // possible first idea to solve the upper note: siccome sto usando solo laplace mi è sufficiente creare un erede di RegressionData. Tuttavia in future implementazioni vorrei poter applicare f-PIRLS anche a problemi con altre condizioni di regolarizzazione (PDE, time varying, etc.) in quel caso probabilmente dovrò costruire appositi figli delle classi derivate.

  std::vector<VectorXr> _solution; //!A Eigen::VectorXr: Stores the system solution.
  std::vector<Real> _dof; //! A vector storing the computed dofs
  std::vector<Real> _J_minima;

    //double scale_parameter = 1; // default value since it is 1 for most of the distributions
    //int scale_parameter_flag = 0; // It tells me if I have to estimate phi (W. do it only for gamma)

  // link and other functions, it can be extended to handle vector of parameters
    virtual Real link(const Real& mu) = 0; // g(.)

    virtual Real link_deriv(const Real& mu) = 0; // g'(.)

    virtual Real inv_link(const Real& theta) = 0; // g^-1(.)

    virtual Real var_function(const Real& mu) = 0; // V(mu)

    virtual void initialize_mu(const VectorXr& y) = 0; // It is different for binary or non-binary outcomes

    void compute_z(); // perform step (1) of f-PIRLS

    void compute_G(); //assemble G matrix

    void compute_W(); // assemble W matrix

    void update_solution(const RegressionData& pseudoData); // perform step (2) of f-PIRLS it is dependent on templates specs

    void compute_mu(const RegressionData& pseudoData); // perform step (3) of f-PIRLS

    bool stopping_criterion(); // it stops PIRLS based on difference between functionals J_k J_k+1 or n_iterations > max_num_iterations

    Real compute_J(); // compute the current value of J

  public:

    FPIRLS_Base(const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& inputData ):mesh_(mesh), inputData_(inputData){}; // Constructor

   void apply(); // perform PIRLS and instanciate the solution in _solution , _dof

   //! A inline member that returns a VectorXr, returns the whole solution_.
   inline std::vector<VectorXr> const & getSolution() const{return _solution;}
   //! A function returning the computed dofs of the model
   inline std::vector<Real> const & getDOF() const{return _dof;}
   //! A function returning the current value of J
   inline std::vector<Real> const & get_J() const{return _J_minima;}


   // FURTHER QUESTIONS - where to implement GCV (lambda can be specified by user or estimated via GCV), where to implement the estimate of phi

   // QUESTION: se avessi un costruttore padre con parametri e non definisco il costruttore figlio (aka, figlio dispone solo di default constr.) e creo un oggetto figlio chiamando il costruttore con parametri ottengo automaticamente quello del padre? In caso negativo devo implementare il costruttore anche per le singole distribuzioni.

};

// Bernoulli
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class FPIRLS_Bernoulli : public FPIRLS_Base <InputHandler, Integrator, ORDER, mydim, ndim> {

  protected:
    inline Real link(const Real& mu){ return log(mu/(1 - mu)); }

    inline Real inv_link(const Real& theta){ return exp(theta)/(1 + exp(theta)); }

    inline Real link_deriv(const Real& mu){ return 1/(mu*(1-mu)); }

    inline Real var_function(const Real& mu){ return(mu*(1-mu)); }

    void initialize_mu(const VectorXr& y) {
      this->mu_.resize(y.size());
      for(UInt i = 0; i < y.size(); i++)
          this->mu_[i] = 1/2 * (y[i] + 1/2); // It is different for binary or non-binary outcomes
      }
  public:

    FPIRLS_Bernoulli(const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& inputData):
      FPIRLS_Base<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData ){};
};


// Poisson
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
    class FPIRLS_Poisson : public FPIRLS_Base <InputHandler, Integrator, ORDER, mydim, ndim> {

    protected:

      inline Real link(const Real& mu){ return log(mu); }

      inline Real link_deriv(const Real& mu){ return 1/mu; }

      inline Real inv_link(const Real& theta){ return exp(theta); }

      inline Real var_function(const Real& mu){ return mu ;}

      void initialize_mu(const VectorXr & y) {
      this->mu_ = y;} // It is different for binary or non-binary outcomes

    public:

    FPIRLS_Poisson(const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& inputData):
      FPIRLS_Base<InputHandler, Integrator, ORDER, mydim, ndim>(mesh, inputData ){};

};

// Exponential
// TO BE DONE

// Gamma
// TO BE DONE

// Probit
// TO BE DONE: how to get normal cdf and its inverse? thay're needed for probit

// Cloglog
// TO BE DONE: W. write his wrong in the old R code, has to be fixed here

#include "FPIRLS_imp.h"

#endif
