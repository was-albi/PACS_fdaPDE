# SIM_GAUSSIAN with right mesh type

#.rs.restartR()
rm(list=ls())
graphics.off()

###############################
## Gam fem fit GAUSSIAN TEST ##
###############################



# LOAD library ----------------------------------------

setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace")
#setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace")
setwd("/Users/giuliopn/PACSworkspace3/PACSworkspace/")
library(fdaPDE)
library(purrr)
#source("GAM_tests/2013_SSR_AllFunctions.R")
#source("GAM_tests/gam_fem_fit.R")
#library(fda)



#### simple mesh 2D (elliptic PDE + covariates + locations different from mesh nodes) ####
{
  
  # Load the mesh
  data(simpleMesh2Ddata)
  
  mesh=create.mesh.2D(nodes=nodes, triangles = triangles, order=2)
  
  plot(mesh)
  
  # Create the FEM basis object
  FEMbasis = create.FEM.basis(mesh)
  
  set.seed(5847947)
  
  # Exact solution
  data_exact=sin(pi*mesh$nodes[,1])
  
  # Plot exact solution
  plot(FEM(data_exact,FEMbasis = FEMbasis))
  
  scale.param = 0.1
  # Perturbed data - locations different from nodes
  data = data_exact + rnorm(n = length(mesh$nodes[,1]), sd = 0.1)
  
  # Set a vector of smoothing coefficients
  lambda = 10^-4
  
  GCVFLAG=FALSE
  GCVMETHODFLAG='Exact'   
  
  # Set PDE parameters
  PDE_parameters_anys = list(K = matrix(c(0.01,0,0,1), nrow = 2), b = c(0,0), c = 0)
  
  # deterministic covariate - Nodes locations
  cov1_nod=sin(pi*mesh$nodes[,1])
  
  #plot covariate 
  image(FEM(cov1_nod,FEMbasis))
  
  cov2_nod=rnorm(mean=0, sd=0.5,n=length(mesh$nodes[,1]))
  W_nod=cbind(cov1_nod,cov2_nod)
  
  
  
  output_CPP <- fdaPDE::gam.fem.fit(observations = data, FEMbasis =FEMbasis, covariates = W_nod, PDE_parameters = PDE_parameters_anys,
                                    lambda = lambda, max.steps=15, fam="gaussian", mu0=NULL, mesh.const=T, scale.param=scale.param)
  
  
  
  output_CPP_classical = smooth.FEM(observations = data, covariates=W_nod,
                          FEMbasis = FEMbasis, lambda = lambda, 
                          PDE_parameters = PDE_parameters_anys,
                          GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)

#image(output_CPP4$fit.FEM)
output_CPP$beta
output_CPP_classical$beta  
  
  
  
}


#### simple mesh 2D (elliptic PDE + covariates + locations different from mesh nodes) ####

rm(list=ls())

# Load the mesh
data(simpleMesh2Ddata)

mesh=create.mesh.2D(nodes=nodes, triangles = triangles, order=2)

plot(mesh)

# Create the FEM basis object
FEMbasis = create.FEM.basis(mesh)

set.seed(5847947)

# Exact solution
data_exact=sin(pi*mesh$nodes[,1])

# Plot exact solution
plot(FEM(data_exact,FEMbasis = FEMbasis))

# Locations different from nodes
xobs=runif(min=-0.5,max=0.5,n=80)
yobs=runif(min=-0.5,max=0.5,n=80)
loc=cbind(xobs,yobs)

# Exact data - locations different from nodes
data_exact = sin(pi*xobs)

# Perturbed data - locations different from nodes
data = data_exact + rnorm(n = length(xobs), sd = 0.1)

# Set a vector of smoothing coefficients
lambda = c(10^-4, 1, 10^4)

lambda = 10^-4


GCVFLAG=FALSE
GCVMETHODFLAG='Exact'   

# Set PDE parameters
PDE_parameters_anys = list(K = matrix(c(0.01,0,0,1), nrow = 2), b = c(0,0), c = 0)

# deterministic covariate - Nodes locations
cov1_nod=sin(pi*mesh$nodes[,1])

#plot covariate 
image(FEM(cov1_nod,FEMbasis))

# Covariates - Locations different from nodes
cov1_nonod=sin(pi*xobs)
cov2_nonod=rnorm(mean=0, sd=0.5,n=length(xobs))
W_nonod=cbind(cov1_nonod,cov2_nonod)



output_CPP <- fdaPDE::gam.fem.fit(observations = data, FEMbasis =FEMbasis, covariates = W_nonod,
                                  lambda = lambda, max.steps=1, fam="gaussian", mu0=NULL, mesh.const=F, scale.param=scale.param)



#image(output_CPP4$fit.FEM)
output_CPP$beta

output_CPP <- fdaPDE::gam.fem.fit(observations = data, locations=loc, FEMbasis =FEMbasis, covariates = W_nonod, PDE_parameters = PDE_parameters_anys,
                                  lambda = lambda, max.steps=15, fam="gaussian", mu0=NULL, mesh.const=F, scale.param=0.1)

output_CPP_classical = smooth.FEM(observations = data, locations=loc, covariates=W_nonod,
                        FEMbasis = FEMbasis, lambda = lambda, 
                        PDE_parameters = PDE_parameters_anys,
                        GCV=GCVFLAG, GCVmethod = GCVMETHODFLAG)

#image(output_CPP4$fit.FEM)
output_CPP$beta
output_CPP_classical$beta
