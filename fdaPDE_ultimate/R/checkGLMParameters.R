checkGLMParameters<-function(max.steps, mu0, mesh.const, method.phi, scale.param, tune, weight, threshold)
{
	  #################### Parameter Check #########################
	# Check max.steps 
	if(!all.equal(max.steps, as.integer(max.steps)) || max.steps < 0 )
		stop("max.steps must be a positive integer.")
	
	if(!is.null(mu0))
  	{
    	if(any(is.na(mu0)))
      		stop("Missing values not admitted in 'mu0'.")

      	if(ncol(mu0) != 1)
    		stop("'mu0' must be a column vector")
  		if(nrow(mu0) < 1)
  			stop("'mu0' must contain at least one element")
  	}

	# check mesh.const
	if (is.null(mesh.const)) 
		stop("'mesh.const' required;  is NULL.")
	
	if(!is.logical(mesh.const))
		stop("'mesh.const' is not logical")
	
	# check method.phi
	if(method.phi != 1 && method.phi!=2)
		stop("'method.phi' must be equal either to 1 or 2")
	
	# check scale.param
	if(!is.null(scale.param)){
		if( !is.numeric(scale.param))
    		stop("'scale.param' must be a real number")
	}
	# check weight
	if(!is.logical(weight))
		stop("'weight' is not logical")

	# check tune
	if (is.null(tune)){ 
		stop("'tune' required;  is NULL.")
	}else if( !is.numeric(tune) || tune < 0){
    	stop("'tune' must be a real positive")
	}
	

	# check threshold
	if (is.null(threshold)){ 
		stop("'threshold' required;  is NULL.")
	}else if( !is.numeric(threshold) || threshold < 0){
    	stop("'threshold' must be a real positive")
	}

}