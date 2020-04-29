/*
 * RegressionOutput.cpp
 *
 *  Created on: Mar 25, 2020
 *      Authors: Colombo - Perin
 */

#define R_VERSION_

#include "fdaPDE.h"
//#include "IO_handler.h"
#include "mesh_objects.h"
#include "mesh.h"
#include "evaluator.h"

extern "C" {
//! This function compute and returns the estimates for Beta and fn in R.
/*!
	This function is then called from R code.

  params:

*/

SEXP get_Regression_Output(SEXP Rmesh, SEXP Rlocations, SEXP Robservations , SEXP Rcovariates, SEXP Rorder ,SEXP RincidenceMatrix,  SEXP RSolCoef, SEXP Rmydim, SEXP Rndim)
{

  MixedFERegressionBase<RegressionData, Integrator, ORDER, mydim, ndim>::computeRegressionSolution(_solution, pseudoData, mesh_, fn_hat, beta_hat);



  		#ifdef R_VERSION_
  		explicit RegressionData(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP Rcovariates,
  								SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF, SEXP RGCVmethod,
  								SEXP Rnrealizations);
  		#endif

}



} // end extern-C
