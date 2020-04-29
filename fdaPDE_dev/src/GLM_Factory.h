#ifndef GLM_FACTORY_
#define GLM_FACTORY_

#include "EXP_Factory.h"

#include <memory>


template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args)
{
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<typename Data_Type,typename Parameter_Type>
class Exp_Family_Factory
{
public:

  static std::unique_ptr<Exp_Family<Data_Type,Parameter_Type>> create_exp_distribution(const std::string &distribution ){

    if(distribution == "Bernoulli") return make_unique<Bernoulli_Dist<Eigen::VectorXi, double>>();
    if(distribution == "Poisson") return make_unique<Poisson_Dist<Eigen::VectorXi, double>>();
    if(distribution == "Normal") return make_unique<Normal_Dist<Eigen::VectorXd, Eigen::Vector2d>>();

  }

};

// TO BE MODIFIED

//! A Factory class: A class for the choice of the cross-validation method to use for the selection of the parameter lambda for each PC.
template<typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
class MixedFEFPCAfactory
{
	public:
	//! A method that takes as parameter a string and builds a pointer to the right object for the cross-validation
	static std::unique_ptr<MixedFEFPCABase<Integrator, ORDER,  mydim,  ndim>> createFPCAsolver(const std::string &validation, const MeshHandler<ORDER,mydim,ndim>& mesh, const FPCAData& fpcaData){
	if(validation=="GCV") return make_unique<MixedFEFPCAGCV<Integrator, ORDER,  mydim, ndim>>(mesh,fpcaData);
	if(validation=="KFold") return make_unique<MixedFEFPCAKFold<Integrator, ORDER,  mydim, ndim>>(mesh,fpcaData);
	if(validation=="NoValidation") return make_unique<MixedFEFPCA<Integrator, ORDER,  mydim, ndim>>(mesh,fpcaData);
	}

};

#endif
