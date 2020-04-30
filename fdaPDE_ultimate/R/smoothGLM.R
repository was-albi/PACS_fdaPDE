#' Spatial generalize linear model with differential regularization
gam.fem.fit= function(locations = NULL, observations, FEMbasis,
	lambda, covariates = NULL, BC = NULL, PDE_parameters=NULL, incidence_matrix = NULL, max.steps=15, mu0=NULL,
	fam=c("binomial", "probit","cloglog", "exponential", "gamma", "poisson", "gaussian"),
	mesh.const=F, method.phi=1,
  scale.param=NULL, GCV = FALSE, nrealizations=100, GCVmethod="Stochastic",
	tune=1.8, weight=F, threshold=0.0004)
{

### 1) Mesh, GVC and FAMILY method set up
### mesh set up
if(class(FEMbasis$mesh) == "mesh.2D"){
    ndim = 2
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.2.5D"){
    ndim = 3
    mydim = 2
  }else if(class(FEMbasis$mesh) == "mesh.3D"){
    ndim = 3
    mydim = 3
  }else{
    stop('Unknown mesh class')
  }
  ### GCVmethod set up
  if(GCVmethod=="Stochastic")
  GCVMETHOD=2
  else if(GCVmethod=="Exact")
  GCVMETHOD=1
  else{
    stop("GCVmethod must be either Stochastic or Exact")
  }
  #### FAMILY set up
  family_admit = c("binomial", "probit", "cloglog", "exponential", "gamma", "poisson", "gaussian")
 if(sum(fam==family_admit)==0 ){
     stop("'family' parameter required.\nCheck if it is one of the following: binomial, probit, cloglog, exponential, gamma, poisson")
}

  # General Check of other parameters
  
  #checkGLMParameters(max.steps, mu0, mesh.const, method.phi, scale.param, tune, weight)
  space_varying=checkSmoothingParameters(locations = locations, observations = observations, FEMbasis= FEMbasis, lambda = lambda, covariates = covariates, incidence_matrix = incidence_matrix, BC = BC, GCV = GCV, PDE_parameters = PDE_parameters, GCVmethod = GCVMETHOD, nrealizations = nrealizations)
  
  checkGLMParameters(max.steps, mu0, mesh.const, method.phi, scale.param, tune, weight, threshold)

  ################## End checking parameters, sizes and conversion #############################
  ## Converting to format internal usage
  if(!is.null(locations))
    locations = as.matrix(locations)
  observations = as.matrix(observations)
  lambda = as.matrix(lambda)
  if(!is.null(covariates))
  covariates = as.matrix(covariates)
  if(!is.null(incidence_matrix))
    incidence_matrix = as.matrix(incidence_matrix)
  if(!is.null(BC))
  {
    BC$BC_indices = as.matrix(BC$BC_indices)
    BC$BC_values = as.matrix(BC$BC_values)
  }

  # if I have PDE non-sv case I need (constant) matrices as parameters
  
  if(!is.null(PDE_parameters) & space_varying==FALSE) 
  {
    PDE_parameters$K = as.matrix(PDE_parameters$K)
    PDE_parameters$b = as.matrix(PDE_parameters$b)
    PDE_parameters$c = as.matrix(PDE_parameters$c)
  }

  checkSmoothingParametersSize(locations, observations, FEMbasis, lambda, covariates, incidence_matrix, BC, GCV, space_varying, PDE_parameters,ndim, mydim)
  
  ################## End checking parameters, sizes and conversion #############################

  if(class(FEMbasis$mesh) == 'mesh.2D' & is.null(PDE_parameters)){

    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.GLM(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
            covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
            BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, FAMILY=fam,
            mu0 = mu0, max.steps=max.steps, mesh.const=mesh.const,
            method.phi=method.phi, scale.param=scale.param, tune=tune, weight=weight, threshold=threshold)

    numnodes = nrow(FEMbasis$mesh$nodes)

    } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==FALSE){
        bigsol = NULL
        print('C++ Code Execution')
        bigsol = CPP_smooth.FEM.GLM.PDE.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
            PDE_parameters = PDE_parameters, 
            covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
            BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, FAMILY=fam,
            mu0 = mu0, max.steps=max.steps, mesh.const=mesh.const,
            method.phi=method.phi, scale.param=scale.param, tune=tune, weight=weight, threshold=threshold )

        numnodes = nrow(FEMbasis$mesh$nodes)
   
   } else if(class(FEMbasis$mesh) == 'mesh.2D' & !is.null(PDE_parameters) & space_varying==TRUE){ 
    
    bigsol = NULL
    print('C++ Code Execution')
    bigsol = CPP_smooth.FEM.GLM.PDE.sv.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
            PDE_parameters = PDE_parameters, 
            covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
            BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, FAMILY=fam,
            mu0 = mu0, max.steps=max.steps, mesh.const=mesh.const,
            method.phi=method.phi, scale.param=scale.param, tune=tune, weight=weight, threshold=threshold )
  
    numnodes = nrow(FEMbasis$mesh$nodes)
  } else if(class(FEMbasis$mesh) == 'mesh.2.5D'){
    
    bigsol = NULL  
    print('C++ Code Execution')
    if(!is.null(locations))
      stop("The option locations!=NULL for manifold domains is currently not implemented")
    bigsol = CPP_smooth.manifold.FEM.GLM.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
            covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
            BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, FAMILY=fam,
            mu0 = mu0, max.steps=max.steps, mesh.const=mesh.const,
            method.phi=method.phi, scale.param=scale.param, tune=tune, weight=weight, threshold=threshold )
    
    numnodes = FEMbasis$mesh$nnodes
    
  }else if(class(FEMbasis$mesh) == 'mesh.3D'){
      
    bigsol = NULL  
    print('C++ Code Execution')
    bigsol = CPP_smooth.volume.FEM.GLM.basis(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
            covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
            BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, FAMILY=fam,
            mu0 = mu0, max.steps=max.steps, mesh.const=mesh.const,
            method.phi=method.phi, scale.param=scale.param, tune=tune, weight=weight, threshold=threshold )
    
    numnodes = FEMbasis$mesh$nnodes
  }
    
  # Estimated functions (smooth fields)
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]# laplacian(f)

  print("DOFs: \n")
  print(bigsol[[2]])
  
  # Make Functional objects object
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)

  reslist = NULL
  beta_hat = bigsol[[4]]
  fn_hat = bigsol[[5]]

  J_minima = bigsol[[3]]

  computedGCV = bigsol[[7]]
  if(GCV == TRUE)
  {
    reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, beta_hat = beta_hat, fn_hat = fn_hat, J_minima = J_minima, GCV = computedGCV)
  }else{
    reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, beta_hat = beta_hat, fn_hat = fn_hat, J_minima = J_minima)
  }
  return(reslist)

  }
