#' Spatial generalize linear model with differential regularization
#'

gam.fem.fit= function(locations = NULL, observations, FEMbasis,
	lambda, covariates = NULL, PDE_parameters=NULL, max.steps=15, mu0=NULL, 
	fam=c("binomial", "probit","cloglog", "exponential", "gamma", "poisson"), 
	mesh.const=T, method.phi=1, 
	psi=NULL, scale.param=NULL, GCV = FALSE, GCVmethod="Stochastic", 
	tune=1.8, weight=F )
{

## ATTENZIONE incidence_matrix e BC non sono variabili ancora implementate, ma apparentemente servono.
##            Nel dubbio le setto inizialmente uguali a NULL, e poi in caso si implementeranno.
incidence_matrix = NULL
BC = NULL

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
  if(fam == "binomial"){
  	FAMILY = 1
    }else if(fam == "probit"){
     FAMILY = 2
     }else if(fam == "cloglog"){
       FAMILY = 3
       }else if(fam == "exponential"){
         FAMILY = 4
         }else if(fam == "gamma"){
           FAMILY = 5
           }else if(fam == "poisson"){
             FAMILY = 6
             }else{
              stop("'fam' required; it can be: binomial, probit, cloglog, exponential, gamma, poisson") 
            }

  # General Check of other parameters
  checkGLMParameters(locations, observations, FEMbasis, lambda, covariates, PDE_parameters, max.steps, mu0, mesh.const, method.phi, psi, scale.param, GCV, tune, weight)

  ################## End checking parameters, sizes and conversion #############################
  ## Converting to format internal usage
  if(!is.null(locations))
  locations = as.matrix(locations)
  observations = as.matrix(observations)
  lambda = as.matrix(lambda)
  if(!is.null(covariates))
  covariates = as.matrix(covariates)

  checkSmoothingParametersSize(locations, observations, FEMbasis, lambda, covariates)

  ################## End checking parameters, sizes and conversion #############################

  if(class(FEMbasis$mesh) == 'mesh.2D' & is.null(PDE_parameters)){	

    bigsol = NULL
    print('C++ Code Execution')
    #ToDo: 2)
    bigsol = CPP_smooth.FEM.GLM(locations=locations, observations=observations, FEMbasis=FEMbasis, lambda=lambda,
      covariates=covariates, incidence_matrix=incidence_matrix, ndim=ndim, mydim=mydim,
      BC=BC, GCV=GCV, GCVMETHOD=GCVMETHOD, nrealizations=nrealizations, FAMILY=FAMILY,
      mu0 = mu0, max.steps=max.steps, mesh.count=mesh.count,
      method.phi=method.phi, psi=psi, scale.param=scale.param, tune=tune, weight=weight )

    numnodes = nrow(FEMbasis$mesh$nodes)
    
    } else{
     stop('At the moment the other cases are not implemented.')
   }
  # Estimated functions (smooth fields)
  f = bigsol[[1]][1:numnodes,]
  g = bigsol[[1]][(numnodes+1):(2*numnodes),]# laplacian(f)

  # Make Functional objects object 
  fit.FEM  = FEM(f, FEMbasis)
  PDEmisfit.FEM = FEM(g, FEMbasis)

  reslist = NULL
  #ToDo: 3) (see 210 in smoothing.R)
  beta = getBetaCoefficients(locations, observations, fit.FEM, covariates, incidence_matrix, ndim, mydim)

  if(GCV == TRUE)
  {
  	#ToDo: 4)
    #seq=getGCV(locations = locations, observations = observations, fit.FEM = fit.FEM, covariates = covariates, incidence_matrix = incidence_matrix, edf = bigsol[[2]], ndim, mydim)
    #reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, beta = beta, edf = bigsol[[2]], stderr = seq$stderr, GCV = seq$GCV)
    }else{
      reslist=list(fit.FEM = fit.FEM, PDEmisfit.FEM = PDEmisfit.FEM, beta = beta)
    }

    return(reslist)

  }