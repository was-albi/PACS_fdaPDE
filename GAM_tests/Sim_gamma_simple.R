# SIM_GAMMA with right mesh type

#.rs.restartR()
rm(list=ls())
graphics.off()

############################
## Gam fem fit GAMMA TEST ##
############################



# LOAD library ----------------------------------------

#setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace")
setwd("/Users/giuliopn/PACSworkspace3/PACSworkspace/")
library(fdaPDE)
library(purrr)
source("GAM_tests/2013_SSR_AllFunctions.R")
source("GAM_tests/gam_fem_fit.R")
library(fda)

####################
# Gam fem fit CPP ##
####################

{
# Build DATA: 2D mesh ---------------------------------
{
  #### square 2D (basic case)
  # setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace/fdaPDE_ultimate/data")
  data(square2Ddata)
  mesh = fdaPDE::create.mesh.2D(nodes=nodes)
  # x11()
  plot(mesh)
  # axis(1)
  # axis(2)
  nnodes = dim(mesh$nodes)[1]
  FEMbasis = fdaPDE::create.FEM.basis(mesh)
}
# Test function ------------------------------------
{
  # Build DATA: 2D mesh ---------------------------------
  {
    #### square 2D (basic case)
    # setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace/fdaPDE_ultimate/data")
    data(square2Ddata)
    mesh = fdaPDE::create.mesh.2D(nodes=nodes)
    # x11()
    plot(mesh)
    # axis(1)
    # axis(2)
    nnodes = dim(mesh$nodes)[1]
    FEMbasis = fdaPDE::create.FEM.basis(mesh)
  }
  # Test function ------------------------------------
  {
    set.seed(5847947)
    # set.seed(42)
    
    a1=runif(1,min=-1.5,max=1.5)
    a2=runif(1,min=-1.5,max=1.5)
    
    z<-function(p)
    {
      
      a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1]) - 2
      
    }
    
    # Exact solution (pointwise at nodes)
    sol_exact=rep(0,dim(mesh$nodes)[1])
    for(i in 1: dim(mesh$nodes)[1])
      sol_exact[i]=z(mesh$nodes[i,])
    
    ran=range(sol_exact) 
    
  }
  
  # Set smoothing parameter ---------------------------
  {
    
    link<-function(x){-1/x}
    inv.link<-function(x){-1/x}
    
    count = 0;
    desmat=matrix(0,nrow=nnodes,ncol=2)
    mu=numeric(length=nnodes)
    scale.param = 1
    
    while(any(mu<=0))
    {
      count <- count+1
      # The seed is set in such a way that if count <10, which is very likely, then
      # the seeds are different for all the simulations (which is highly desirable!)
      set.seed(42 + count)
      desmat[,1]=rbeta(nnodes,shape1=1.5,shape2=2)
      desmat[,2]=rbeta(nnodes,shape1=3,shape2=2)+1
      beta1=-2/5
      beta2=3/10
      betas_truth = c(beta1,beta2)
      param=sol_exact+beta1*desmat[,1]+beta2*desmat[,2]
      mu<-inv.link(param)
      #if(all(mu>=0)) seeds.check[sim] <- seeds[sim] + count
    }
    
    response <- rgamma(nnodes, shape=mu/scale.param, scale=scale.param)
  }
  
  GCVFLAG=T
  GCVMETHODFLAG='Exact'
  lambda= 10^-4
  
  #mu_guessed <- rep(1,nnodes) + response
  
  
  output_CPP <- fdaPDE::gam.fem.fit(observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                                    lambda = lambda, max.steps=10, fam="gamma", mu0=NULL, mesh.const=T, scale.param=scale.param)
  
  betas_hat_CPP = output_CPP$beta_hat
  #Best result for lambda: 10^-2, 0.05
  
  # loro
  #image(FEM(sol_exact, FEMbasis))
  # nostra
  #image(FEM(output_CPP$fit.FEM$coeff,FEMbasis))
  
  ########################
  # Gam fem fit Wilhelm ##
  ########################
  
  # Set up mesh object for Whilelm function
  {
    # parm: nodes(or triangular vertix), edge=NULL, 
    basisobj = create.FEM.basis(mesh$nodes, e = NULL, mesh$triangles, mesh$order)
    
    #  set up a dummy FEM functional data object
    simfd = fd(numeric(FEMbasis$nbasis),basisobj)
  }
  # set up observation class for Whilelm function
  observed.data = cbind(1:nnodes,as.numeric(response))
  np = nnodes # it must be initialize in Whilelm function
  # set lambda
  lambda_W = lambda
  
  logit <- function(x){qlogis(x)}
  inv.logit <- function(x){plogis(x)}
  
  setwd("/home/alb/Scrivania/PACS/Git_Folder/")
  
  Whilelm.fit = gam.fem.fit( data = observed.data ,desmat ,fdobj = simfd ,lambda = lambda_W, max.steps=10, 
                             fam="gamma", mesh.const=T, method.phi=1, psi=NULL, 
                             scale.param=scale.param, GCV.score=FALSE, tune=1.8, weight=F )
  
  Whilelm.fit$beta
  #image(FEM(Whilelm.fit$felsplobj$coefs,FEMbasis))
  
  betas_hat_W = Whilelm.fit$beta
  
  # Result ---------------------------------
  format(betas_hat_W, digits = 15)
  format(betas_hat_CPP,digits = 15)
  
  #  Percentage error
  (betas_hat_W - betas_truth)/betas_truth * 100 
  (betas_hat_CPP - betas_truth)/betas_truth * 100 
  
  
  plot(Whilelm.fit$felsplobj$coefs)
  plot(output_CPP$fn_hat)
  plot(sol_exact)
  
  
  
  library(ggplot2)
  
  f_data = data.frame(Whilelm.fit$felsplobj$coefs,output_CPP$fn_hat)   
  names(f_data) <- c("Whilelm_fn_hat","CPP_fn_hat")
  
  ggplot(f_data, aes(y=Whilelm_fn_hat, x = CPP_fn_hat)) + geom_point() + 
    geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5) +
    ggtitle("Comparison solution coefficients") +
    xlab("CPP f coefficients") +
    ylab("R f coefficients") +
    theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  element_name_CPP <- "apply_return_laplacian_coefs_.txt"
  
  setwd("/home/alb/Scrivania/PACS/Git_Folder/debugging_output/CPP/")
  my_lapl <- read.table(element_name_CPP, header = FALSE, sep = ",", dec = ".", row.names = NULL)
  
  
  lapl_data = data.frame(Whilelm.fit$laplacefd$coefs,my_lapl)   
  names(lapl_data) <- c("Whilelm_lapl","CPP_lapl")
  
  ggplot(lapl_data, aes(y=Whilelm_lapl, x = CPP_lapl)) + geom_point() + 
    geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5) +
    ggtitle("Comparison laplacian coefficients") +
    xlab("CPP Df coefficients") +
    ylab("R Df coefficients") +
    theme(plot.title = element_text(size = 14, hjust = 0.5))
  
  
  image(FEM(Whilelm.fit$felsplobj$coefs,FEMbasis)) # whil
  box()
  title(main = "R solution", font.main = 4)
  image(FEM(output_CPP$fit.FEM$coeff,FEMbasis)) # ours
  title(main ="CPP solution", font.main = 4)
  
  image(FEM(sol_exact, FEMbasis))
  
  
}

GCVFLAG=T
GCVMETHODFLAG='Exact'
lambda= 10^-4



output_CPP <- fdaPDE::gam.fem.fit(observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                                  lambda = lambda, max.steps=10, fam="gamma", mu0=NULL, mesh.const=T, scale.param=scale.param)


#### simple mesh 2D (elliptic PDE + covariates + locations different from mesh nodes) ####
#.rs.restartR()
rm(list=ls())
graphics.off()

#setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace")
setwd("/Users/giuliopn/PACSworkspace3/PACSworkspace/")
library(fdaPDE)
library(purrr)


link <- function(x){-1/x}
inv.link  <- function(x){-1/x}

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
sol_exact = sin(pi*xobs)

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

# Covariates - Locations different from nodes
cov1_nonod=sin(pi*xobs)
cov2_nonod=rnorm(mean=5, sd=0.5,n=length(xobs))
W_nonod=cbind(cov1_nonod,cov2_nonod)

Whilelm.fit = gam.fem.fit( data = observed.data ,desmat ,fdobj = simfd ,lambda = lambda_W, max.steps=10, 
                           fam="gamma", mesh.const=T, method.phi=1, psi=NULL, 
                           scale.param=scale.param, GCV.score=FALSE, tune=1.8, weight=F )

beta1=-2/5
beta2=-3/10
betas_truth = c(beta1,beta2)
param=sol_exact+beta1*W_nonod[,1]+beta2*W_nonod[,2]
#param = sol_exact
mu<-inv.link(param)

scale.param = 1
response <- rgamma(length(xobs), shape=mu/scale.param, scale=scale.param)


# Result ---------------------------------
format(betas_hat_W, digits = 15)
format(betas_hat_CPP,digits = 15)

output_CPP <- fdaPDE::gam.fem.fit(observations = as.numeric(response), locations = loc, FEMbasis =FEMbasis, covariates = W_nonod,
                                  lambda = lambda, max.steps=15, fam="gamma", mu0=NULL, mesh.const=F, scale.param=scale.param)


plot(Whilelm.fit$felsplobj$coefs)
plot(output_CPP$fn_hat)
plot(sol_exact)



library(ggplot2)

f_data = data.frame(Whilelm.fit$felsplobj$coefs,output_CPP$fn_hat)   
names(f_data) <- c("Whilelm_fn_hat","CPP_fn_hat")

ggplot(f_data, aes(y=Whilelm_fn_hat, x = CPP_fn_hat)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5) +
  ggtitle("Comparison solution coefficients") +
  xlab("CPP f coefficients") +
  ylab("R f coefficients") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))

element_name_CPP <- "apply_return_laplacian_coefs_.txt"

setwd("/home/alb/Scrivania/PACS/Git_Folder/debugging_output/CPP/")
my_lapl <- read.table(element_name_CPP, header = FALSE, sep = ",", dec = ".", row.names = NULL)


lapl_data = data.frame(Whilelm.fit$laplacefd$coefs,my_lapl)   
names(lapl_data) <- c("Whilelm_lapl","CPP_lapl")

ggplot(lapl_data, aes(y=Whilelm_lapl, x = CPP_lapl)) + geom_point() + 
  geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5) +
  ggtitle("Comparison laplacian coefficients") +
  xlab("CPP Df coefficients") +
  ylab("R Df coefficients") +
  theme(plot.title = element_text(size = 14, hjust = 0.5))


image(FEM(Whilelm.fit$felsplobj$coefs,FEMbasis)) # whil
box()
title(main = "R solution", font.main = 4)
image(FEM(output_CPP$fit.FEM$coeff,FEMbasis)) # ours
title(main ="CPP solution", font.main = 4)

image(FEM(sol_exact, FEMbasis))


}






#### simple mesh 2D (elliptic PDE + covariates + locations different from mesh nodes) ####
#.rs.restartR()
rm(list=ls())
graphics.off()

setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace")
library(fdaPDE)
library(purrr)


link <- function(x){-1/x}
inv.link  <- function(x){-1/x}

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
sol_exact = sin(pi*xobs)

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

# Covariates - Locations different from nodes
cov1_nonod=sin(pi*xobs)
cov2_nonod=rnorm(mean=5, sd=0.5,n=length(xobs))
W_nonod=cbind(cov1_nonod,cov2_nonod)

beta1=-2/5
beta2=-3/10
betas_truth = c(beta1,beta2)
param=sol_exact+beta1*W_nonod[,1]+beta2*W_nonod[,2]
#param = sol_exact
mu<-inv.link(param)

scale.param = 1
response <- rgamma(nnodes, shape=mu/scale.param, scale=scale.param)



output_CPP <- fdaPDE::gam.fem.fit(observations = as.numeric(response), FEMbasis =FEMbasis, covariates = W_nonod,
                                  lambda = lambda, max.steps=1, fam="gamma", mu0=NULL, mesh.const=F, scale.param=scale.param)


#image(output_CPP4$fit.FEM)
output_CPP$beta
