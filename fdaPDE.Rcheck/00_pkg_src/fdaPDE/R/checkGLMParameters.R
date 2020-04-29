checkGLMParameters<-function(locations = NULL, observations, FEMbasis,
	lambda, covariates = NULL, PDE_parameters=NULL, max.steps=15, mu0=NULL, 
	fam=c("binomial", "probit","cloglog", "exponential", "gamma", "poisson"), 
	mesh.const=T, method.phi=1, 
	psi=NULL, scale.param=NULL, GCV = FALSE, GCVmethod="Stochastic", 
	tune=1.8, weight=F )
{
	#################### Parameter Check #########################
	# check locations
	if(!is.null(locations))
	{
		if(any(is.na(locations)))
		stop("Missing values not admitted in 'locations'.")
		if(any(is.na(observations)))
		stop("Missing values not admitted in 'observations' when 'locations' are specified.")
	}
	# check observations
	if (is.null(observations))
	stop("observations required;  is NULL.")
	# check FEMbasis
	if (is.null(FEMbasis))
	stop("FEMbasis required;  is NULL.")
	if(class(FEMbasis)!= "FEMbasis")
	stop("'FEMbasis' is not class 'FEMbasis'")

	if(class(FEMbasis$mesh)!='MESH2D' & class(FEMbasis$mesh) != "MESH.2.5D" & class(FEMbasis$mesh) != "MESH.3D")
	stop('Unknown mesh class')

	if((class(FEMbasis$mesh) == "MESH.2.5D" || class(FEMbasis$mesh) == "MESH.3D") & !is.null(PDE_parameters) )
	stop('For mesh classes different from MESH2D, anysotropic regularization is not yet implemented. 
		Use Laplacian regularization instead')

	# check lambda
	if (is.null(lambda))
	stop("lambda required;  is NULL.")

	# GLM model implement only the laplacian case. In the future it could be implemented.
	if(!is.null(PDE_parameters))
	{
		stop("At the moment GLM model support only Laplacian case. 'PDE_parameters' must be NULL")
	}
	# Check max.steps 
	if(!all.equal(max.steps, as.integer(max.steps)) || max.steps < 0 )
	stop("max.steps must be a positive integer.")

	# check mesh.const
	if (is.null(mesh.const)) 
	stop("'mesh.const' required;  is NULL.")
	if(!is.logical(mesh.const))
	stop("'mesh.const' is not logical")
	# check method.phi
	if(method.phi != 1 || method.phi!=2)
	stop("'method.phi' must be equal either to 1 or 2")

	# check GCV
	if (is.null(GCV)) 
	stop("GCV required;  is NULL.")
	if(!is.logical(GCV))
	stop("'GCV' is not logical")

}

checkGLMParametersSize<-function(locations = NULL, observations, FEMbasis, lambda, covariates = NULL)
{
	if(ncol(observations) != 1)
	stop("'observations' must be a column vector")
	if(nrow(observations) < 1)
	stop("'observations' must contain at least one element")

	if(is.null(locations))
	{
		if(class(FEMbasis$mesh) == "MESH2D"){
			if(nrow(observations) > nrow(FEMbasis$mesh$nodes))
			stop("Size of 'observations' is larger then the size of 'nodes' in the mesh")
			}else if(class(FEMbasis$mesh) == "MESH.2.5D" || class(FEMbasis$mesh) == "MESH.3D"){
				if(nrow(observations) > FEMbasis$mesh$nnodes)
				stop("Size of 'observations' is larger then the size of 'nodes' in the mesh")
			}
		}
		if(ncol(lambda) != 1)
		stop("'lambda' must be a column vector")
		if(nrow(lambda) < 1)
		stop("'lambda' must contain at least one element")
		if(!is.null(covariates))
		{
			if(nrow(covariates) != nrow(observations))
			stop("'covariates' and 'observations' have incompatible size;")
		}

}