#ifndef __REGRESSIONESTIMATES_HPP__
#define __REGRESSIONESTIMATES_HPP__

#include "fdaPDE.h"
#include "mesh.h"
#include "regressionData.h"
#include "evaluator.h"

#include "printer.h"


/*! A class for computation of beta coefficient and function estimation
*/

template<typename InputHandler, UInt ORDER, UInt mydim, UInt ndim>
class RegressionEstimates{

	VectorXr _observations;
	MeshHandler<ORDER,mydim,ndim> _mesh;
	std::vector<Point> _points;// location points
	MatrixXr _CovMatrix;
	VectorXr _WeightsMatrix;
	//MatrixXr _IncidenceMatrix;

	// coefficient expressing the solution w.r.t basis of the Finite Element space
	std::vector<VectorXr> _solution_coeff;

	// evaluation of the function in the observation points
	std::vector<VectorXr> fn_hat_;
	std::vector<VectorXr> beta_hat_;

	#ifdef R_VERSION_
		void setPoints(SEXP Rlocations);
		void setObservations(SEXP Robservations);
		void setSolutionCoeff(SEXP Rsolution);
		void setCovariates(SEXP Rcovariates);
		void setWeightsMatrix(SEXP RweightsMatrix);
		// void setIncidenceMatrix(SEXP RincidenceMatrix);
	#endif

public:

	// R constructr
	#ifdef R_VERSION_
	explicit RegressionEstimates(SEXP Rlocations, SEXP Robservations, MeshHandler<ORDER,mydim,ndim>& mesh, SEXP Rcovariates, SEXP RweightsMatrix, SEXP Rsolution);
	#endif

	// c++ constructr
	explicit RegressionEstimates(const InputHandler& regressionData, const MeshHandler<ORDER,mydim,ndim>& mesh, const std::vector<VectorXr>& solution);

	explicit RegressionEstimates(const InputHandler& regressionData, const MeshHandler<ORDER,mydim,ndim>& mesh, const VectorXr& solution);


	//! A method computing the result of the regression (beta_hat, fn_hat) given the solution of the FEM problem
	void computeEstimates();

	// get methods
	inline std::vector<VectorXr> const& getFunctionEst()const{return fn_hat_;};
	inline std::vector<VectorXr> const& getBetaEst()const{return beta_hat_;};


};


#include "regressionEstimates_imp.h"

#endif
