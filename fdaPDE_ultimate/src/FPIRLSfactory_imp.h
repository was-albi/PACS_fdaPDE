#ifndef __FPIRLSAFACTORY_HPP__
#define __FPIRLSAFACTORY_HPP__

#include"FPIRLSfactory.h"

/*
template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
std::unique_ptr<FPIRLS_Base<InputHandler, Integrator, ORDER,  mydim,  ndim>> FPIRLSfactory<InputHandler, Integrator, ORDER, mydim, ndim>::createFPIRLSsolver(const std::string &family, const MeshHandler<ORDER,mydim,ndim>& mesh, const InputHandler& inputData, VectorXr& mu0){

if(mu0.size() == 0){
  VectorXr y = inputData.getObservations();
  if(family == "binomial"){ //binary outcomes
    mu0 = VectorXr::Zero(y.size());
    for(UInt i = 0; i < y.size(); i++){
        mu0[i] = 0.5 * (y[i] + 0.5);
    }//end for-i
  }else{ // not-binary outcomes
    mu0 = y;
  }
}//end if

if(family=="binomial"){
    return make_unique<FPIRLS_Bernoulli<InputHandler, Integrator, ORDER,  mydim, ndim>>(mesh,inputData,mu0);
}else if(family=="poisson"){
    return make_unique<FPIRLS_Poisson<InputHandler, Integrator, ORDER,  mydim, ndim>>(mesh,inputData,mu0);
}else if(family=="gamma"){
    return make_unique<FPIRLS_Gamma<InputHandler, Integrator, ORDER,  mydim, ndim>>(mesh,inputData,mu0);
}

}


#ifdef R_VERSION_

template <typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
	VectorXr FPIRLSfactory<InputHandler, Integrator, ORDER, mydim,ndim>::convert_mu(SEXP Rmu0){

    VectorXr mu0;
    UInt n_obs_ = Rf_length(Rmu0);
  	mu0.resize(n_obs_);

  	UInt count = 0;
  	for(UInt i=0;i<n_obs_;++i)
  	   mu0[i] = REAL(Rmu0)[i];

    return mu0;
  }


#endif
*/
#endif
