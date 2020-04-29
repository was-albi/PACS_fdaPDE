
// Same constructor except for WeightsMatrix
RegressionData::RegressionData(std::vector<Point>& locations, VectorXr& observations, UInt order, std::vector<Real> lambda, MatrixXr& covariates, MatrixXr& WeightsMatrix, MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices, std::vector<Real>& bc_values, bool DOF):
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




RegressionDataGAM::RegressionDataGAM(std::vector<Point>& locations, VectorXr& observations, UInt order,
                std::vector<Real> lambda, UInt max_num_iterations, Real threshold, MatrixXr& covariates,
                MatrixXi& incidenceMatrix, std::vector<UInt>& bc_indices,
                std::vector<Real>& bc_values, bool DOF):
		 RegressionData(locations, observations, order, lambda, covariates, incidenceMatrix, bc_indices, bc_values, DOF), max_num_iterations_(max_num_iterations), threshold_(threshold)
{;}


#ifdef R_VERSION_

RegressionDataGAM::RegressionDataGAM(SEXP Rlocations, SEXP Robservations, SEXP Rorder, SEXP Rlambda, SEXP Rmax_num_iteration, SEXP Rthreshold, SEXP Rcovariates, SEXP RincidenceMatrix, SEXP RBCIndices, SEXP RBCValues, SEXP DOF,SEXP RGCVmethod, SEXP Rnrealizations):
	RegressionData(Rlocations, Robservations, Rorder, Rlambda, Rcovariates, RincidenceMatrix,
					 			   RBCIndices, RBCValues, DOF,RGCVmethod, Rnrealizations)
{
	max_num_iterations_ = INTEGER(Rmax_num_iteration)[0];
	threshold_ =  REAL(Rthreshold)[0];
}

#endif

RegressionData RegressionDataGAM::generatePseudodata(VectorXr& input_pseudo_Obs, VectorXr& W )const{
	// Copy Input params:
	std::vector<Point> input_location = this->getLocations();
	UInt input_order = this->getOrder();
	std::vector<Real> input_lambda = this->getLambda();
	MatrixXi input_incidenceMatrix = this->getIncidenceMatrix();
	std::vector<UInt> input_bc_indices = this->getDirichletIndices();
	std::vector<Real> input_bc_values = this->getDirichletValues();
	bool input_dof = this->computeDOF();
	MatrixXr input_Cov = this -> getCovariates();

  const UInt W_size = W.size();
	MatrixXr input_WeightsMatrix = MatrixXr::Zero(W_size,W_size);
	for(UInt i = 0; i<W_size; i++)
		input_WeightsMatrix(i,i) = W(i);

	RegressionData pseudoData ( input_location, input_pseudo_Obs,  input_order,  input_lambda, input_Cov, input_WeightsMatrix, input_incidenceMatrix, input_bc_indices, input_bc_values , input_dof);

	return pseudoData;
}
