 template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
	void computeRegressionSolution(const std::vector<VectorXr>& solution,const InputHandler& regressionData ,std::vector<VectorXr>&  fn_hat ,std::vector<VectorXr>&  beta_hat){

    std::vector<Point> points = regressionData.getLocations();
    VectorXr solution_coeff;
    const UInt Points_Length = points.size();
    std::vector<bool> isinside(Points_Length);
    MatrixXr WeightsMatrix = regressionData_.getWeightsMatrix();
    MatrixXr CovMatrix = regressionData_.getCovariates();
    Evaluator<ORDER,2,2> f_evaluator(mesh_);

    solution_coeff.resize(_solution[0].size()/2);

    for(UInt j = 0; j < regressionData.getLambda().size(); j++ ){

      fn_hat[j].resize(Points_Length);

      for(UInt i = 0; i < _solution[0].size()/2; i++){ // take only the value of f and not of laplacian(f)
         solution_coeff(i) = _solution[j](i);
      } //end for-i

      f_evaluator.eval(points, solution_coeff, true, fn_hat[j], isinside);

      // solve weighted least squares
      beta_hat[j] = (WeightsMatrix*CovMatrix).bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(WeightsMatrix*(Observations - fn_hat[j])) ;

    } //end of for-j
}// end method





 template<typename InputHandler, typename Integrator, UInt ORDER, UInt mydim, UInt ndim>
 void MixedFERegressionBase<InputHandler,Integrator,ORDER,mydim, ndim>::setH()
 {

 	UInt nlocations = regressionData_.getNumberofObservations();


 	MatrixXr W(this->regressionData_.getCovariates());

 	if(regressionData_.isLocationsByNodes())
 	{
 		MatrixXr W_reduced(regressionData_.getNumberofObservations(), W.cols());
 		for (auto i=0; i<nlocations;++i)
 		{
 			auto index_i = regressionData_.getObservationsIndices()[i];
 			for (auto j=0; j<W.cols();++j)
 			{
 				W_reduced(i,j) = W(index_i,j);
 			}
 		}
 		W = W_reduced;
 	}

	MatrixXr WeightsMatrix(this->regressionData_.getWeightsMatrix());

 if(WeightsMatrix.rows()==0){ // if regression is not weighted

	 MatrixXr WTW(W.transpose()*W);
   H_=W*WTW.ldlt().solve(W.transpose()); // using cholesky LDLT decomposition for computing hat matrix

 }else{ // if regression IS weighted

	 MatrixXr WTW(W.transpose()*WeightsMatrix*W);
	 H_=W*WTW.ldlt().solve(W.transpose()*WeightsMatrix); // using cholesky LDLT decomposition for computing hat matrix

 }

}
