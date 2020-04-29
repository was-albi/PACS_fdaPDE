
#define R_VERSION_

#include "fdaPDE.h"
#include "regressionData.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "finite_element.h"
#include "matrix_assembler.h"
#include "FPCAData.h"
#include "FPCAObject.h"
#include "solverdefinitions.h"
//#include <chrono>

#include "mixedFEFPCA.h"
#include "mixedFERegression.h"
#include "mixedFEFPCAfactory.h"


// The skeleton actually create the solver and call the PIRLS algo, then build up the SEXP variable to be returned in R
// REMARK: it is not embedded in a extern C scope, since it has not to be called by R even if return SEXP, indeed it is called by gam_fem_fit (see below)

template<typename InputHandler,typename Integrator,UInt ORDER, UInt mydim, UInt ndim>
SEXP GAM_skeleton(InputHandler &GAMData, SEXP Rmesh, std::string family)
{

  MeshHandler<ORDER, mydim, ndim> mesh(Rmesh);

// Factory:
  // std::unique_ptr<FPRILS_Base<InputHandler, Integrator, ORDER, mydim, ndim>> fpirls = FPIRLSfactory<Integrator, ORDER, mydim, ndim>::createPIRLSsolver(family, mesh, GAMData);

  FPIRLS_Bernoulli<InputHandler, Integrator, ORDER, mydim, ndim> fpirls(mesh, GAMData);

  fpirls->apply();

  const std::vector<VectorXr>& solution = fpirls.getSolution();
	const std::vector<Real>& dof = fpirls.getDOF();
  const std::vector<Real> J_value = fpirls.get_J();

// COMPOSIZIONE SEXP result FOR RETURN

	//Copy result in R memory
	SEXP result = NILSXP;
	result = PROTECT(Rf_allocVector(VECSXP, 3));
	SET_VECTOR_ELT(result, 0, Rf_allocMatrix(REALSXP, solution[0].size(), solution.size()));
	SET_VECTOR_ELT(result, 1, Rf_allocVector(REALSXP, solution.size()));
  SET_VECTOR_ELT(result, 2, Rf_allocVector(REALSXP, J_value.size()))
	Real *rans = REAL(VECTOR_ELT(result, 0));
	for(UInt j = 0; j < solution.size(); j++)
	{
		for(UInt i = 0; i < solution[0].size(); i++)
			rans[i + solution[0].size()*j] = solution[j][i];
	}

	Real *rans2 = REAL(VECTOR_ELT(result, 1));
	for(UInt i = 0; i < solution.size(); i++)
	{
		rans2[i] = dof[i];
	}

  Real *rans3 = REAL(VECTOR_ELT(result, 2));
  for(UInt i = 0; i < J_value.size(); i++)
	{
		rans2[i] = J_value[i];
	}

	UNPROTECT(1);
	return(result);
}




extern "C" {

//! This function manages the various options for Spatial Regression, Sangalli et al version
/*!
	This function is then called from R code.
	\param Robservations an R-vector containing the values of the observations.
	\param Rdesmat an R-matrix containing the design matrix for the regression.
	\param Rmesh an R-object containg the output mesh from Trilibrary
	\param Rorder an R-integer containing the order of the approximating basis.
	\param Rlambda an R-double containing the penalization term of the empirical evidence respect to the prior one.
	\param Rcovariates an R-matrix of covariates for the regression model
	\param RincidenceMatrix an R-matrix containing the incidence matrix defining the regions for the smooth regression with areal data
	\param RBCIndices an R-integer containing the indexes of the nodes the user want to apply a Dirichlet Condition,
			the other are automatically considered in Neumann Condition.
	\param RBCValues an R-double containing the value to impose for the Dirichlet condition, on the indexes specified in RBCIndices
	\param DOF an R boolean indicating whether dofs of the model have to be computed or not
	\param RGCVmethod an R-integer indicating the method to use to compute the dofs when DOF is TRUE, can be either 1 (exact) or 2 (stochastic)
	\param Rnrealizations the number of random points used in the stochastic computation of the dofs
	\return R-vector containg the coefficients of the solution
*/

// SEXP regression_Laplace(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
	//				SEXP Rlambda, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
	//				SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations)

  SEXP gam_fem_fit(SEXP Rlocations, SEXP Robservations, SEXP Rmesh, SEXP Rorder,SEXP Rmydim, SEXP Rndim,
  					SEXP Rlambda, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues,
  					SEXP DOF, SEXP RGCVmethod, SEXP Rnrealizations , SEXP Rfamily, SEXP Rmax_num_iteration, SEXP Rtreshold, SEXP tune, SEXP mu0, SEXP phi, SEXP scale_param )
{
    //Set input data: can I re-use this fdapde class? FPCA actually implement it's own class for the datastructure
	RegressionDataGAM regressionData(Rlocations, Robservations, Rorder, Rlambda, Rmax_num_iteration, Rtreshold, Rcovariates, RincidenceMatrix, RBCIndices, RBCValues, DOF, RGCVmethod, Rnrealizations);

 // mydim e ndim sono entrambi 2 nei nostri casi (mesh 2D)
 // ndim si riferisce allo spazio in cui la mesh è embedded, mydim è la dimensione della mesh (in caso Rimannian manifold si ha 2d immerso in 3d)
	UInt mydim=INTEGER(Rmydim)[0];
	UInt ndim=INTEGER(Rndim)[0];

  std::string family = CHAR(STRING_ELT(Rfamily,0));

    if(regressionData.getOrder()==1 && mydim==2 && ndim==2)
    	return(GAM_skeleton<RegressionData,IntegratorTriangleP2, 1, 2, 2>(regressionData, Rmesh, family));
    else if(regressionData.getOrder()==2 && mydim==2 && ndim==2)
		return(GAM_skeleton<RegressionData,IntegratorTriangleP4, 2, 2, 2>(regressionData, Rmesh, family));
    //else if(regressionData.getOrder()==1 && mydim==2 && ndim==3)
    //"NOT IMPLEMENTED YET"
    //return(regression_skeleton<RegressionData,IntegratorTriangleP2, 1, 2, 3>(regressionData, Rmesh));
   // else if(regressionData.getOrder()==2 && mydim==2 && ndim==3)
		//"NOT IMPLEMENTED YET"
    //return(regression_skeleton<RegressionData,IntegratorTriangleP4, 2, 2, 3>(regressionData, Rmesh));
	//else if(regressionData.getOrder()==1 && mydim==3 && ndim==3)
    // "NOT IMPLEMENTED YET"
    //return(regression_skeleton<RegressionData,IntegratorTetrahedronP2, 1, 3, 3>(regressionData, Rmesh));
    return(NILSXP);
}

} // end extern C scope
