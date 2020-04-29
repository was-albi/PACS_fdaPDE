CPP_smooth.FEM.GLM<-function(locations, observations, FEMbasis, lambda,
                                  covariates=NULL, incidence_matrix=NULL, ndim, mydim,
                                  BC=NULL, GCV, GCVMETHOD = 2, nrealizations=100, mu0=NULL, max.steps=15, FAMILY, 
                                  mesh.const, method.phi=1, psi=NULL, scale.param=NULL, tune=1.8, weight=F)
{
  # Indexes in C++ starts from 0, in R from 1, opportune transformation
  FEMbasis$mesh$triangles = FEMbasis$mesh$triangles - 1
  FEMbasis$mesh$edges = FEMbasis$mesh$edges - 1
  FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] = FEMbasis$mesh$neighbors[FEMbasis$mesh$neighbors != -1] - 1
  if(is.null(covariates))
  {
    covariates<-matrix(nrow = 0, ncol = 1)
  }

  if(is.null(locations))
  {
    locations<-matrix(nrow = 0, ncol = 2)
  }
  
  if(is.null(incidence_matrix))
  {
    incidence_matrix<-matrix(nrow = 0, ncol = 1)
  }
  
  if(is.null(BC$BC_indices))
  {
    BC$BC_indices<-vector(length=0)
  }else
  {
    BC$BC_indices<-as.vector(BC$BC_indices)
  }
  
  if(is.null(BC$BC_values))
  {
    BC$BC_values<-vector(length=0)
  }else
  {
    BC$BC_values<-as.vector(BC$BC_values)
  }

  if(is.null(mu0))
  {
    mu0<-matrix(nrow = 0, ncol = 1)
  }
  if(is.null(psi))
  {
    psi<-matrix(nrow = 0, ncol = 1)
  }
  if(is.null(scale.param))
  {
    scale.param<-matrix(nrow = 0, ncol = 1)
  }


  ## Set proper type for correct C++ reading
  locations <- as.matrix(locations)
  storage.mode(locations) <- "double"
  storage.mode(FEMbasis$mesh$nodes) <- "double"
  storage.mode(FEMbasis$mesh$triangles) <- "integer"
  storage.mode(FEMbasis$mesh$edges) <- "integer"
  storage.mode(FEMbasis$mesh$neighbors) <- "integer"
  storage.mode(FEMbasis$order) <- "integer"
  covariates <- as.matrix(covariates)
  storage.mode(covariates) <- "double"
  incidence_matrix <- as.matrix(incidence_matrix)
  storage.mode(incidence_matrix) <- "integer"
  storage.mode(ndim) <- "integer"
  storage.mode(mydim) <- "integer"
  storage.mode(lambda) <- "double"
  storage.mode(BC$BC_indices) <- "integer"
  storage.mode(BC$BC_values) <-"double"
  
  GCV <- as.integer(GCV)
  storage.mode(GCV) <-"integer"
  
  storage.mode(nrealizations) <- "integer"
  storage.mode(GCVMETHOD) <- "integer"

  storage.mode(FAMILY) <- "integer"
  storage.mode(max.steps) <- "integer"
  storage.mode(tune) <- "double"
  storage.mode(mu0) <- "double"
  storage.mode(psi) <- "double"
  storage.mode(scale.param) <- "double"
  
  # Treshold is not chosen by the users in W. code, but a fixed value is used
  threshold = 10^(-4)
  storage.mode(treshold) <- "double"
  
  # Test for FAMILY as string
  FAMILY = "bernoulli"
  storage.mode(FAMILY) <- "character"
  
  
  
  ## Call C++ function

  bigsol <- .Call("gam_fem_fit", locations, observations, FEMbasis$mesh, FEMbasis$order,
                  mydim, ndim, lambda, covariates, incidence_matrix, BC$BC_indices, BC$BC_values,
                  GCV, GCVMETHOD, nrealizations, FAMILY, max.steps, tune, mu0, psi, scale.param, PACKAGE = "fdaPDE")
  # GCV became DOF variable in C++, the motivation is if GCV = TRUE we need to compute DoF otherwise no. 
  return(bigsol)

}