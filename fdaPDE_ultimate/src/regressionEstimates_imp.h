#ifndef __REGRESSIONESTIMATES_IMP_HPP__
#define __REGRESSIONESTIMATES_IMP_HPP__

#include "regressionEstimates.h"


template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
	RegressionEstimates<InputHandler, ORDER, mydim, ndim>::RegressionEstimates(const InputHandler& regressionData,
										const MeshHandler<ORDER,mydim,ndim>& mesh, const std::vector<VectorXr>& solution):_mesh(mesh)
{

	// class members initialization
	 _observations = regressionData.getObservations();
	 _points = regressionData.getLocations();

	 _CovMatrix = regressionData.getCovariates();
	 _WeightsMatrix = regressionData.getWeightsMatrix();

	 //initialize _solution_coeff
	 VectorXr temp_solution = VectorXr::Zero(solution[0].size()/2);
	 for(UInt j = 0; j < regressionData.getLambda().size(); j++ ){
	 	for(UInt i = 0; i < solution[0].size()/2; i++){ // take only the value of f and not of laplacian(f)
				temp_solution(i) = solution[j](i);
	 	} //end for-i
	 	_solution_coeff.push_back(temp_solution);

	}//end for-j

}


#ifdef R_VERSION_
template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
	RegressionEstimates<InputHandler, ORDER, mydim, ndim>::RegressionEstimates(SEXP Rlocations, SEXP Robservations, MeshHandler<ORDER,mydim,ndim>& mesh,
										SEXP Rcovariates, SEXP RweightsMatrix, SEXP Rsolution):_mesh(mesh)
{
	setPoints(Rlocations);
	//setIncidenceMatrix(RincidenceMatrix);
	setObservations(Robservations);
	setCovariates(Rcovariates);
	setWeightsMatrix(RweightsMatrix);
	setSolutionCoeff(Rsolution);

}

#endif



//! A method computing the result of the regression (beta_hat, fn_hat) given the solution of the FEM problemi
template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
	void RegressionEstimates<InputHandler, ORDER, mydim, ndim>::computeEstimates(){

	 const UInt Points_Length = _points.size();
	 const UInt length = _solution_coeff.size();
	 std::vector<bool> isinside(Points_Length);
	 Evaluator<ORDER, mydim, ndim> f_evaluator(_mesh);

	 VectorXr temp;

	 for(UInt j = 0; j < length; j++ ){

		if(Points_Length >0 ){ // Location points are specified, different from mesh nodes

		 	temp.resize(Points_Length);
		 	f_evaluator.eval(_points, _solution_coeff[j], true, temp, isinside);
		 	fn_hat_.push_back(temp);

		 }else{// Locations coicide with mesh nodes
		 	fn_hat_.push_back( _solution_coeff[j]);
		 }


		 if(_CovMatrix.rows()>0){
			 //MatrixXr Cov = _CovMatrix;
			 MatrixXr Cov = _CovMatrix;
			 VectorXr _obs_centered = _observations - fn_hat_[j];

				 MatrixXr WTW = Cov.transpose()*_WeightsMatrix.asDiagonal()*Cov;
				 //VectorXr rightterm = Cov.transpose()*_WeightsMatrix*_obs_centered;
				 MatrixXr desmatprod = WTW.colPivHouseholderQr().solve(Cov.transpose()*_WeightsMatrix.asDiagonal());

				 //VectorXr tmp = WTW.colPivHouseholderQr().solve(rightterm);
				 //VectorXr tmp = desmatprod*_obs_centered;

				 VectorXr tmp = VectorXr::Zero(desmatprod.rows());
			/*	 Real current_product_value = 0;
				 for(UInt i = 0; i < desmatprod.rows();i++){
         		for(UInt j=0; j< desmatprod.cols();j++){
							current_product_value = desmatprod(i,j)*_obs_centered(j);
							if(abs(current_product_value) < 0.00000001)
								current_product_value = 0;
							tmp(i) += current_product_value;
						}
				 }*/
				 tmp = desmatprod*_obs_centered;
				 beta_hat_.push_back(tmp);
	 	 }

/*
		 std::string saving_filename = "regressionestimates_betahat_";
		 saving_filename = saving_filename + ".txt";
		 printer::saveVectorXr(saving_filename,beta_hat_[0]);
*/

	 } //end of for-j
}// end method


#ifdef R_VERSION_

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void RegressionEstimates<InputHandler, ORDER, mydim, ndim>::setPoints(SEXP Rlocations)
{
	UInt n = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[0];
	if(n>0){
		int pointsDim = INTEGER(Rf_getAttrib(Rlocations, R_DimSymbol))[1];

	  if (pointsDim == 2){
			for(auto i=0; i<n; ++i)
			{
				_points.emplace_back(REAL(Rlocations)[i+ n*0],REAL(Rlocations)[i+ n*1]);
			}
		}else{ //ndim == 3
			for(auto i=0; i<n; ++i)
			{
				_points.emplace_back(REAL(Rlocations)[i+ n*0],REAL(Rlocations)[i+ n*1],REAL(Rlocations)[i+ n*2]);
			}
		}
	}
}

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void RegressionEstimates<InputHandler, ORDER, mydim, ndim>::setObservations(SEXP Robservations)
{
	UInt n_obs_ = Rf_length(Robservations);
	_observations.resize(n_obs_);
	//observations_indices_.reserve(n_obs_);

	UInt count = 0;
	if(_points.size() == 0 ) //&& nRegions_ == 0
	{
		//locations_by_nodes_ = true;
		for(UInt i=0;i<n_obs_;++i)
		{
			if(!ISNA(REAL(Robservations)[i]))
			{
				_observations[count] = REAL(Robservations)[i];
				count++;
				//observations_indices_.push_back(i);
			}
		}
		_observations.conservativeResize(count, Eigen::NoChange);
	}
	else // locations_.size() > 0 NOR nRegions_ > 0
	{
		//locations_by_nodes_ = false;
		for(UInt i=0;i<n_obs_;++i)
			_observations[i] = REAL(Robservations)[i];
	}
}

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void RegressionEstimates<InputHandler, ORDER, mydim, ndim>::setCovariates(SEXP Rcovariates)
{
	UInt n_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[0];
	UInt p_ = INTEGER(Rf_getAttrib(Rcovariates, R_DimSymbol))[1];

	_CovMatrix.resize(n_, p_);

	for(UInt i=0; i<n_; ++i)
	{
		for(UInt j=0; j<p_ ; ++j)
			_CovMatrix(i,j)=REAL(Rcovariates)[i+ n_*j];
	}
}

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void RegressionEstimates<InputHandler, ORDER, mydim, ndim>::setWeightsMatrix(SEXP RweightsMatrix){
	UInt Psize = INTEGER(Rf_getAttrib(RweightsMatrix, R_DimSymbol))[0];

	_WeightsMatrix = VectorXr::Zero(Psize);

	for(UInt i=0; i<Psize; ++i)
			_WeightsMatrix(i)=REAL(RweightsMatrix)[i];

}

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void RegressionEstimates<InputHandler, ORDER, mydim, ndim>::setSolutionCoeff(SEXP Rsolution){
	// Rsolution collect vector solution by column
	UInt nRows = INTEGER(Rf_getAttrib(Rsolution, R_DimSymbol))[0];
	UInt nCols = INTEGER(Rf_getAttrib(Rsolution, R_DimSymbol))[1];

	_solution_coeff.resize(nCols);

	for(UInt j = 0; j< nCols ; j++){
		_solution_coeff[j].resize(nRows);
		for(UInt i = 0; i< nRows ; i++)
			_solution_coeff[j](i) = REAL(Rsolution)[i + j*nRows ];
	}
}
#endif

/*
template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
void RegressionEstimates<InputHandler, ORDER, mydim, ndim>::setIncidenceMatrix(SEXP RincidenceMatrix)
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
*/
#endif
