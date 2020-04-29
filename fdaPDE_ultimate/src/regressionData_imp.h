#ifndef __REGRESSIONDATA_IMP_HPP__
#define __REGRESSIONDATA_IMP_HPP__

RegressionData::RegressionData(std::vector<Point>& locations, VectorXr& observations, UInt order, std::vector<Real> lambda, MatrixXr& covariates, MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF):
					locations_(locations), observations_(observations), covariates_(covariates), incidenceMatrix_(incidenceMatrix),
					order_(order), lambda_(lambda), bc_values_(bc_values), bc_indices_(bc_indices), DOF_(DOF)
{
	nRegions_ = incidenceMatrix_.rows();
	if(locations_.size()==0 && nRegions_==0)
	{
		locations_by_nodes_ = true;
		for(int i = 0; i<observations_.size();++i) observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}
}

// Same constructor except for WeightsMatrix
RegressionData::RegressionData(std::vector<Point>& locations, VectorXr& observations, UInt order, std::vector<Real> lambda, MatrixXr& covariates, VectorXr& WeightsMatrix, MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF):
					locations_(locations), observations_(observations), covariates_(covariates), WeightsMatrix_(WeightsMatrix), incidenceMatrix_(incidenceMatrix),
					order_(order), lambda_(lambda), bc_values_(bc_values), bc_indices_(bc_indices), DOF_(DOF)
{
	nRegions_ = incidenceMatrix_.rows();
	if(locations_.size()==0 && nRegions_==0)
	{
			locations_by_nodes_ = true;
			for(int i = 0; i<observations_.size();++i) observations_indices_.push_back(i);
	}
	else
	{
		locations_by_nodes_ = false;
	}
}


RegressionDataElliptic::RegressionDataElliptic(std::vector<Point>& locations, VectorXr& observations, UInt order,
												std::vector<Real> lambda, Eigen::Matrix<Real,2,2>& K,
												Eigen::Matrix<Real,2,1>& beta, Real c, MatrixXr& covariates,
												MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,
												std::vector<Real>& bc_values, bool DOF):
		 RegressionData(locations, observations, order, lambda, covariates, incidenceMatrix, bc_indices, bc_values, DOF), K_(K), beta_(beta), c_(c)
{;}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(std::vector<Point>& locations,
									VectorXr& observations, UInt order, std::vector<Real> lambda,
									const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > >& K,
									const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > >& beta,
									const std::vector<Real>& c, const std::vector<Real>& u,
									MatrixXr& covariates, MatrixXi& incidenceMatrix,
									std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF):
		RegressionData(locations, observations, order, lambda, covariates, incidenceMatrix, bc_indices, bc_values, DOF), K_(K), beta_(beta), c_(c), u_(u)
{;}


#ifdef R_VERSION_
RegressionData::RegressionData(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP Rcovariates,
							SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF, SEXP RGCVmethod,
							SEXP Rnrealizations)
{
	setLocations(Rlocations);
	setIncidenceMatrix(RincidenceMatrix);
	setObservations(Robservations);
	setCovariates(Rcovariates);
	setNrealizations(Rnrealizations);

	GCVmethod_ = INTEGER(RGCVmethod)[0];

	order_ =  INTEGER(Rorder)[0];
	DOF_ = INTEGER(DOF)[0];
	UInt length_indexes = Rf_length(RBCIndices);

	bc_indices_.assign(INTEGER(RBCIndices), INTEGER(RBCIndices) +  length_indexes);

	bc_values_.assign(REAL(RBCValues),REAL(RBCValues) + Rf_length(RBCIndices));

    UInt length_lambda = Rf_length(Rlambda);
    for (UInt i = 0; i<length_lambda; ++i)  lambda_.push_back(REAL(Rlambda)[i]);

}

RegressionDataElliptic::RegressionDataElliptic(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta,
				 SEXP Rc, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations):
	RegressionData(Rlocations, Robservations, Rorder, Rlambda, Rcovariates, RincidenceMatrix,
					 			   RBCIndices, RBCValues, DOF,RGCVmethod, Rnrealizations)
{
	K_.resize(2, 2);
	for(auto i=0; i<2; ++i)
	{
		for(auto j=0; j<2 ; ++j)
		{
			K_(i,j)=REAL(RK)[i+ 2*j];
		}
	}

	beta_.resize(2);
	for(auto i=0; i<2 ; ++i)
	{
		beta_(i)=REAL(Rbeta)[i];
	}

	c_ =  REAL(Rc)[0];
}

RegressionDataEllipticSpaceVarying::RegressionDataEllipticSpaceVarying(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta,
				 SEXP Rc, SEXP Ru, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations):
					 RegressionData(Rlocations, Robservations, Rorder, Rlambda, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, DOF,RGCVmethod, Rnrealizations),
					 K_(RK), beta_(Rbeta), c_(Rc), u_(Ru)
{;}

#endif

void RegressionDataEllipticSpaceVarying::print(std::ostream & out) const
{
	for (auto i=0;i<18;i++)
		out<<K_(i);
	for (auto i=0;i<18;i++)
		out<<beta_(i);
	for (auto i=0;i<18;i++)
		out<<c_(i);
}

// create GAM constructors with right template inheritances
// Laplace
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(std::vector<Point>& locations, VectorXr& observations, UInt order,
                std::vector<Real> lambda, UInt max_num_iterations, Real threshold, MatrixXr& covariates,
                MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,
                std::vector<Real>& bc_values, bool DOF, Real tune):
		 RegressionData(locations, observations, order, lambda, covariates, incidenceMatrix, bc_indices, bc_values, DOF), max_num_iterations_(max_num_iterations), threshold_(threshold), tune_param(tune)
{;}

// PDE
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(std::vector<Point>& locations, VectorXr& observations, UInt order,
								std::vector<Real> lambda, Eigen::Matrix<Real,2,2>& K,
								Eigen::Matrix<Real,2,1>& beta, Real c, UInt max_num_iterations, Real threshold, MatrixXr& covariates,
								MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,
								std::vector<Real>& bc_values, bool DOF, Real tune):
		 RegressionDataEllipticSpaceVarying(locations, observations, order, lambda, K, beta, c, covariates, incidenceMatrix, bc_indices, bc_values, DOF), max_num_iterations_(max_num_iterations), threshold_(threshold), tune_param(tune)
{;}

// PDE SpaceVarying
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(std::vector<Point>& locations, VectorXr& observations,
											UInt order, std::vector<Real> lambda,
											const std::vector<Eigen::Matrix<Real,2,2>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,2> > >& K,
											const std::vector<Eigen::Matrix<Real,2,1>, Eigen::aligned_allocator<Eigen::Matrix<Real,2,1> > >& beta,
											const std::vector<Real>& c, const std::vector<Real>& u, UInt max_num_iterations, Real threshold,
											MatrixXr& covariates, MatrixXi& incidenceMatrix,
											std::vector<UInt>& bc_indices, std::vector<Real>& bc_values,
											bool DOF, Real tune):
		 RegressionDataEllipticSpaceVarying(locations, observations, order, lambda, K, beta, c, u, covariates, incidenceMatrix, bc_indices, bc_values, DOF), max_num_iterations_(max_num_iterations), threshold_(threshold), tune_param(tune)
{;}


#ifdef R_VERSION_

// Laplace
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF,SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rtune):
	RegressionData(Rlocations, Robservations, Rorder, Rlambda, Rcovariates, RincidenceMatrix,
					 			   RBCIndices, RBCValues, DOF,RGCVmethod, Rnrealizations)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	tune_param = REAL(Rtune)[0];
	//initialObservations_(observations_);
//  setInitialObservations(Robservations);
}

// PDE
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF,SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rtune):
	RegressionDataElliptic(Rlocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc , Rcovariates, RincidenceMatrix,
					 			   RBCIndices, RBCValues, DOF,RGCVmethod, Rnrealizations)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	tune_param = REAL(Rtune)[0];
	//initialObservations_(observations_);
//	setInitialObservations(Robservations);
}

// PDE SpaceVarying
template<typename RegressionHandler>
RegressionDataGAM<RegressionHandler>::RegressionDataGAM(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP RK, SEXP Rbeta, SEXP Rc, SEXP Ru, SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF,SEXP RGCVmethod, SEXP Rnrealizations, SEXP Rtune):
	RegressionDataEllipticSpaceVarying(Rlocations, Robservations, Rorder, Rlambda, RK, Rbeta, Rc, Ru, Rcovariates, RincidenceMatrix,
					 			   RBCIndices, RBCValues, DOF,RGCVmethod, Rnrealizations)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
	tune_param = REAL(Rtune)[0];
	//initialObservations_(observations_);
//	setInitialObservations(Robservations);
}


#endif

template<typename RegressionHandler>
void RegressionDataGAM<RegressionHandler>::updatePseudodata(VectorXr& z_,  VectorXr& P){

	this-> observations_ = z_;
	this-> WeightsMatrix_ = P;

}

#ifdef R_VERSION_
/*
template<typename RegressionHandler>
void RegressionDataGAM<RegressionHandler>::setInitialObservations(SEXP Robservations)
{
	UInt n_obs_ = Rf_length(Robservations);
	initialObservations_.resize(n_obs_);
	initial_observations_indeces_.reserve(n_obs_);

	UInt count = 0;
	if(this->locations_.size() == 0 && this->nRegions_ == 0)
	{
		this->locations_by_nodes_ = true;
		for(auto i=0;i<n_obs_;++i)
		{
			if(!ISNA(REAL(Robservations)[i]))
			{
				initial_observations_indeces_[count] = REAL(Robservations)[i];
				count++;
				initial_observations_indeces_.push_back(i);
			}
		}
		initialObservations_.conservativeResize(count, Eigen::NoChange);
	}
	else // locations_.size() > 0 NOR nRegions_ > 0
	{
		this->locations_by_nodes_ = false;
		for(auto i=0;i<n_obs_;++i)
		{
			initialObservations_[i] = REAL(Robservations)[i];
		}
	}

}
*/
void RegressionData::setObservations(SEXP Robservations)
{
	UInt n_obs_ = Rf_length(Robservations);
	observations_.resize(n_obs_);
	observations_indices_.reserve(n_obs_);

	UInt count = 0;
	if(locations_.size() == 0 && nRegions_ == 0)
	{
		locations_by_nodes_ = true;
		for(auto i=0;i<n_obs_;++i)
		{
			if(!ISNA(REAL(Robservations)[i]))
			{
				observations_[count] = REAL(Robservations)[i];
				count++;
				observations_indices_.push_back(i);
			}
		}
		observations_.conservativeResize(count, Eigen::NoChange);
	}
	else // locations_.size() > 0 NOR nRegions_ > 0
	{
		locations_by_nodes_ = false;
		for(auto i=0;i<n_obs_;++i)
		{
			observations_[i] = REAL(Robservations)[i];
		}
	}

	//std::cout<<"Observations #"<<observations_.size()<<std::endl<<observations_<<std::endl;
	//for(auto i=0;i<observations_indices_.size();++i)	std::cout<<observations_indices_[i]<<std::endl;
}

void RegressionData::setCovariates(SEXP Rcovariates)
{
	n_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[0];
	p_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[1];

	covariates_.resize(n_, p_);

	for(auto i=0; i<n_; ++i)
	{
		for(auto j=0; j<p_ ; ++j)
		{
			covariates_(i,j)=REAL(Rcovariates)[i+ n_*j];
		}
	}
}

void RegressionData::setLocations(SEXP Rlocations)
{
	n_ = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	if(n_>0){
		int ndim = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[1];

	  if (ndim == 2){
			for(auto i=0; i<n_; ++i)
			{
				locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1]);
			}
		}else{ //ndim == 3
			for(auto i=0; i<n_; ++i)
			{
				locations_.emplace_back(REAL(Rlocations)[i+ n_*0],REAL(Rlocations)[i+ n_*1],REAL(Rlocations)[i+ n_*2]);
			}
		}
	}
}

void RegressionData::setNrealizations(SEXP Rnrealizations) {
	nrealizations_ = INTEGER(Rnrealizations)[0];
}

void RegressionData::setIncidenceMatrix(SEXP RincidenceMatrix)
{
	nRegions_ = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[0];
	UInt p = INTEGER(Rf_getAttrib(RincidenceMatrix, R_DimSymbol))[1];

	incidenceMatrix_.resize(nRegions_, p);

	for(auto i=0; i<nRegions_; ++i)
	{
		for(auto j=0; j<p; ++j)
		{
			incidenceMatrix_(i,j) = INTEGER(RincidenceMatrix)[i+nRegions_*j];
		}
	}
}

#endif

void RegressionData::printObservations(std::ostream & out) const
{

	for(auto i=0;i<observations_.size(); i++)
	{
		out<<i<<"\t"<<observations_(i)<<std::endl;
	}
}

void RegressionData::printCovariates(std::ostream & out) const
{

	for(auto i=0;i<covariates_.rows(); i++)
	{
		for(auto j=0; j<covariates_.cols(); j++)
		{
			out<<covariates_(i,j)<<"\t";
		}
		out<<std::endl;
	}
}

void RegressionData::printLocations(std::ostream & out) const
{

	for(std::vector<Point>::size_type i=0;i<locations_.size(); i++)
	{
		locations_[i].print(out);
		//std::cout<<std::endl;
	}
}

void RegressionData::printIncidenceMatrix(std::ostream & out) const
{
	for (auto i=0; i<incidenceMatrix_.rows(); i++)
	{
		for (auto j=0; j<incidenceMatrix_.cols(); j++)
		{
			out << incidenceMatrix_(i,j) << "\t";
		}
		out << std::endl;
	}
}

#endif
