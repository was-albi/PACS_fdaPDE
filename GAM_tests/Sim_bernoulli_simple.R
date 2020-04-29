# SIM_BERNOULLI with right mesh type

#.rs.restartR()
rm(list=ls())
graphics.off()

################################  
## Gam fem fit BERNOULLI TEST ##
################################


# LOAD library ----------------------------------------

setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace")
#setwd("/Users/giuliopn/PACSworkspace3/PACSworkspace/")
library(fdaPDE)
library(purrr)
source("GAM_tests/2013_SSR_AllFunctions.R")
source("GAM_tests/gam_fem_fit.R")
source("GAM_tests/gam_fem.R")
library(fda)


logit <- function(x){qlogis(x)}
inv.logit <- function(x){plogis(x)}


####################
# Gam fem fit CPP ##
####################

# Build DATA: 2D mesh ---------------------------------
{
  #### square 2D (basic case)
  # setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace/fdaPDE_ultimate/data")
  data(square2Ddata)
  mesh = fdaPDE::create.mesh.2D(nodes=nodes)
  # x11()
  plot(mesh, lwd=3, cex = 1.9)
  # axis(1)
  # axis(2)
  nnodes = dim(mesh$nodes)[1]
  FEMbasis = fdaPDE::create.FEM.basis(mesh)
}
# Test function ------------------------------------

set.seed(5847947)

####################
# Gam fem fit CPP ##
####################
{
  # Build DATA: 2D mesh ---------------------------------
  {
    #### square 2D (basic case)
    setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace/fdaPDE_ultimate/data")
    #setwd("/Users/giuliopn/PACSworkspace3/PACSworkspace/fdaPDE_ultimate/data")
    data(square2Ddata)
    mesh = fdaPDE::create.mesh.2D(nodes=nodes)
    # x11()
    plot(mesh, lwd=3, cex = 1.9)
    # axis(1)
    # axis(2)
    nnodes = dim(mesh$nodes)[1]
    FEMbasis = fdaPDE::create.FEM.basis(mesh)
  }
  # Test function ------------------------------------
  {
    set.seed(5847947)
    
    a1=runif(1,min=-1.5,max=1.5)
    a2=runif(1,min=-1.5,max=1.5)
    
    z<-function(p)
    {
      a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1])
      
    }
    
    # Exact solution (pointwise at nodes)
    sol_exact=rep(0,dim(mesh$nodes)[1])
    for(i in 1: dim(mesh$nodes)[1])
      sol_exact[i]=z(mesh$nodes[i,])
    
    ran=range(sol_exact) 
    
  }
  
  # Set smoothing parameter ---------------------------
  {
    
    link<-function(x){logit(x)}
    inv.link<-function(x){1/(1+exp(-x))}
    
    
    desmat=matrix(0,nrow=nnodes,ncol=2)
    
    desmat[,1]=rbeta(nnodes,shape1=1.5,shape2=2)  # sampling covariates from beta distr.
    desmat[,2]=rbeta(nnodes,shape1=3,shape2=2)+1  # sampling covariates from beta distr.
    beta1=-2/5
    beta2=3/10
    betas_truth = c(beta1,beta2)
    param=sol_exact+beta1*desmat[,1]+beta2*desmat[,2]
    #param = sol_exact
    mu<-inv.link(param)
    response <- rbernoulli(nnodes,p = mu)
  }
  
  #lambda= seq(10^-6,10^-3,by=5*10^-6)
  lambda= 10^-2
  lambda = c(10^-2, 1, 10^2)
  
  GCVFLAG=T
  GCVMETHODFLAG='Exact'
  
  #mu_guessed <- rep(0.5,nnodes)
  
  output_CPP <- fdaPDE::gam.fem.fit(observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat, GCV=T, GCVmethod = 'Exact',
                                    lambda = lambda, max.steps=15, fam="binomial", mu0=NULL, mesh.const=T, scale.param=NULL)
  
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
  
  covariates<-matrix(nrow = 0, ncol = nnodes)
  
  setwd("/home/alb/Scrivania/PACS/Git_Folder/")
  #setwd("/Users/giuliopn/PACSworkspace3/PACSworkspace/GAM_tests/")
  file_name = "debugging_output/R/covariates.txt"
  write.table(format(desmat,scientific = TRUE), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
  
  file_name = "debugging_output/R/link_response.txt"
  tmp_resp <- link(as.numeric(response))
  write.table(format(tmp_resp,scientific = TRUE), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
  
  
  Whilelm.fit = gam.fem.fit( data = observed.data ,desmat=desmat ,fdobj = simfd ,lambda = lambda_W, max.steps=15, 
                             fam="binomial", mesh.const=T, method.phi=1, psi=NULL, 
                             scale.param=NULL, GCV.score=FALSE, tune=1.8, weight=F )
  
  Whilelm.fit = gam.fem( data = observed.data ,desmat=desmat ,fdobj = simfd ,log.lambda.seq = lambda_W, max.steps=15, 
                         fam="binomial", mesh.const=T, method.phi=1, psi=NULL, show=F,
                         scale.param=NULL, tune=1.8 )
  
  
  Whilelm.fit$beta
  # image(FEM(Whilelm.fit$felsplobj$coefs,FEMbasis))
  
  betas_hat_W = Whilelm.fit$beta
  
  # Result --------------------------------- 
  format(betas_hat_W, digits = 15)
  format(betas_hat_CPP,digits = 15)
  #  Percentage error
  (betas_hat_W - betas_truth)/betas_truth * 100 
  (betas_hat_CPP - betas_truth)/betas_truth * 100 
  
  plot(Whilelm.fit$felsplobj$coefs)
  plot(output_CPP$fn_hat)
  
  format(Whilelm.fit$felsplobj$coefs, digits = 15)
  
  format(output_CPP$fn_hat,digits=15)
  
  format(abs(output_CPP$fn_hat-Whilelm.fit$felsplobj$coefs),digits=15)
  
  
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
  
  
  b1_R = 1.045180822994545
  b2_R = 0.886476957854119
  
  b1_CPP = 1.045180822994481
  b2_CPP = 0.886476957853976
  
  b1_R - b1_CPP
  b2_R - b2_CPP
  



output_CPP <- fdaPDE::gam.fem.fit(observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat,
                          lambda = lambda, max.steps=1, fam="binomial", mu0=NULL, mesh.const=T, scale.param=NULL)

#### simple mesh 2D (elliptic PDE + covariates + locations different from mesh nodes) ####
#.rs.restartR()
rm(list=ls())
graphics.off()

setwd("/Users/giuliopn/PACSworkspace3/PACSworkspace/")
library(fdaPDE)
library(purrr)
source("GAM_tests/2013_SSR_AllFunctions.R")
source("GAM_tests/gam_fem_fit.R")
library(fda)


link <- function(x){qlogis(x)}
inv.link <- function(x){plogis(x)}

# Load the mesh
data(simpleMesh2Ddata)

mesh=create.mesh.2D(nodes=nodes, triangles = triangles, order=2)

plot(mesh)

# Create the FEM basis object
FEMbasis = create.FEM.basis(mesh)

set.seed(5847947)

covariates<-matrix(nrow = 0, ncol = nnodes)
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


# Result --------------------------------- 
 format(betas_hat_W, digits = 15)
 format(betas_hat_CPP,digits = 15)
  #  Percentage error
(betas_hat_W - betas_truth)/betas_truth * 100 
(betas_hat_CPP - betas_truth)/betas_truth * 100 
 
   plot(Whilelm.fit$felsplobj$coefs)
   plot(output_CPP$fn_hat)
   
   format(Whilelm.fit$felsplobj$coefs, digits = 15)
   
   format(output_CPP$fn_hat,digits=15)
  
   format(abs(output_CPP$fn_hat-Whilelm.fit$felsplobj$coefs),digits=15)
   
    
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
  
#plot covariate 
image(FEM(cov1_nod,FEMbasis))

# Covariates - Locations different from nodes
cov1_nonod=sin(pi*xobs)
cov2_nonod=rnorm(mean=0, sd=0.5,n=length(xobs))
W_nonod=cbind(cov1_nonod,cov2_nonod)

beta1=-2/5
beta2=3/10
betas_truth = c(beta1,beta2)
param=sol_exact+beta1*W_nonod[,1]+beta2*W_nonod[,2]
#param = sol_exact
mu<-inv.link(param)
response <- rbernoulli(nnodes,p = mu)


output_CPP <- fdaPDE::gam.fem.fit(observations = as.numeric(response), FEMbasis =FEMbasis, covariates = W_nonod,
                                  lambda = lambda, max.steps=1, fam="binomial", mu0=NULL, mesh.const=F, scale.param=NULL)


#image(output_CPP4$fit.FEM)
output_CPP$beta

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


b1_R = 1.045180822994545
b2_R = 0.886476957854119

b1_CPP = 1.045180822994481
b2_CPP = 0.886476957853976

b1_R - b1_CPP
b2_R - b2_CPP

}


#### simple mesh 2D (elliptic PDE + covariates + locations different from mesh nodes) ####
#.rs.restartR()
rm(list=ls())
graphics.off()

setwd("/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace")
library(fdaPDE)
library(purrr)
source("GAM_tests/2013_SSR_AllFunctions.R")
source("GAM_tests/gam_fem_fit.R")
library(fda)


link <- function(x){qlogis(x)}
inv.link <- function(x){plogis(x)}

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
cov2_nonod=rnorm(mean=0, sd=0.5,n=length(xobs))
W_nonod=cbind(cov1_nonod,cov2_nonod)

beta1=-2/5
beta2=3/10
betas_truth = c(beta1,beta2)
param=sol_exact+beta1*W_nonod[,1]+beta2*W_nonod[,2]
#param = sol_exact
mu<-inv.link(param)
response <- rbernoulli(nnodes,p = mu)


output_CPP <- fdaPDE::gam.fem.fit(observations = as.numeric(response), FEMbasis =FEMbasis, covariates = W_nonod,
                                  lambda = lambda, max.steps=1, fam="binomial", mu0=NULL, mesh.const=F, scale.param=NULL)


#image(output_CPP4$fit.FEM)
output_CPP$beta









