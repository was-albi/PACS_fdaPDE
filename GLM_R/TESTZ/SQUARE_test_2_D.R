# TESTING INTERFACE SQUARE 2D

#.rs.restartR()
rm(list=ls())
graphics.off()


# PATH ------------------------------------------------
gam_fem_path = "/home/alb/Scrivania/PACS/Codes_GAM_FEM/"
# gam_fem_path = "/Users/giuliopn/PACSworkspace4/PACS_fdaPDE/"
gamfemfitfolder = "/home/alb/Scrivania/PACS/Git_Folder/"
# gamfemfitfolder = (DOVE VOGLIAMO SALVARE GLI AUTPUT DEI PRINT IN GAMFEMFIT)
mesh_path = "/home/alb/Scrivania/PACS/Git_Folder/PACSworkspace/fdaPDE_ultimate/data"

# LOAD library ----------------------------------------

setwd(gam_fem_path)
library(fdaPDE)
library(purrr)
source("2013_SSR_AllFunctions.R")
source("gam_fem_fit.R")
source("gam_fem.R")
library(fda)

# FAMILY CHOICE ----------------------------------------

FAMILY = "binomial" # "cloglog", "probit", "poisson", "gamma", "exponential" , ( "gaussian", "inv_gaussian" )

{
if(FAMILY == "binomial"){
  logit <- function(x){qlogis(x)}
  inv.logit <- function(x){plogis(x)}
  link = logit
  inv.link = inv.logit
}

if(FAMILY == "cloglog"){
  l<-make.link("cloglog")
  link<-l$linkfun
  inv.link<-l$linkinv
}

if(FAMILY == "probit"){
  l<-make.link("probit")
  link<-l$linkfun
  inv.link<-l$linkinv
}

if(FAMILY == "poisson"){
  l<-make.link("log")
  link<-l$linkfun
  inv.link<-l$linkinv
}


if(FAMILY == "exponential"){
  link<-function(x){-1/x}
  inv.link<-link 
}


if(FAMILY == "gamma"){
  link<-function(x){-1/x}
  inv.link<-link
}

if(FAMILY == "gaussian"){
  link<-function(x){x}
  inv.link<-link
}

if(FAMILY == "invgaussian"){
  link<-function(x){-1/(2*x*x)}
  inv.link<-function(x){ sqrt(-1/(2*x))}
}

}

# BETA ------------------------------------------------

beta1=-2/5
beta2=3/10
betas_truth = c(beta1,beta2)

# lambda ---------------------------------------------

lambda = c(0.01,1,100)
GCVFLAG=T
GCVmethod='Exact'

# scale param -----------------------------------------

scale.param = 1

# Confront with whilhelm? -----------------------------

gamfemfitto = TRUE;

# covariates ------------------------------------------

covariates_flag = TRUE;

# location != nodes -----------------------------------

location_flag = FALSE;

# PDE params ------------------------------------------
# NOT IMPLEMENTED YET 
is_PDE = FALSE;
is_spacevarying = FALSE;

# -------------------- CODE ---------------------------

# mesh reading:
setwd(mesh_path)

data(square2Ddata)

mesh = fdaPDE::create.mesh.2D(nodes=nodes)
# x11()
plot(mesh, lwd=3, cex = 1.9)
# axis(1)
# axis(2)
nnodes = dim(mesh$nodes)[1]
FEMbasis = fdaPDE::create.FEM.basis(mesh)


# locations -----------------------------------------------
set.seed(42)
# locations
if(location_flag){
  # Locations different from nodes
  xobs=runif(min=0,max=0.1,n=80)
  yobs=runif(min=0,max=0.1,n=80)
  loc=cbind(xobs,yobs)
}else{
  loc = mesh$nodes
}


# 2D random field (function f) ------------------------------
a1=runif(1,min=-1.5,max=1.5)
a2=runif(1,min=-1.5,max=1.5)

z<-function(p)
{
  a1*sin(2*pi*p[1])*cos(2*pi*p[2])+a2*sin(3*pi*p[1]) -2
  
}

# exact solution
sol_exact=rep(0,length(loc[,1]))
for(i in 1:length(loc[,1])){
  sol_exact[i] <- z(loc[i,])
}


ran=range(sol_exact) 


# covariates ---------------------------------------------------------

if(covariates_flag){
  
  desmat=matrix(0,nrow=length(loc[,1]),ncol=2)
  
  desmat[,1]=rbeta(length(loc[,1]),shape1=1.5,shape2=2)  # sampling covariates from beta distr.
  desmat[,2]=rbeta(length(loc[,1]),shape1=3,shape2=2)+1  # sampling covariates from beta distr.
  param=sol_exact+beta1*desmat[,1]+beta2*desmat[,2]
}else{
  param = sol_exact
  desmat = NULL
}


mu<-inv.link(param)

  # sampling response:
{
  if(FAMILY == "binomial" || FAMILY == "cloglog" || FAMILY == "probit"){
    response <- rbernoulli(length(loc[,1]),p = mu)
  }
  
  if(FAMILY == "gamma" ){
    response <- rgamma(length(loc[,1]), shape=mu/scale.param, scale=scale.param)
  }

  if(FAMILY == "exponential"){
    response <- rgamma(length(loc[,1]), alpha = 1, scale = mu)
  }
  
  if(FAMILY == "poisson"){
    response <- rpois(length(loc[,1]), lambda = mu)
  }
  
  if(FAMILY=="gaussian"){
    response <- rnorm(length(loc[,1]),mu,scale.param)
  }

  if(FAMILY == "invgaussian"){
    response <- rinvgauss(length(loc[,1]), mu, lambda=scale.param)
  }
  }
  

output_CPP <- fdaPDE::gam.fem.fit(location = loc, observations = as.numeric(response), FEMbasis =FEMbasis, covariates = desmat, GCV=GCVFLAG, GCVmethod = GCVmethod,
                                  lambda = lambda, max.steps=15, fam=FAMILY, mu0=NULL, mesh.const=T, scale.param=NULL)
  
betas_hat_CPP = output_CPP$beta_hat  




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
observed.data = cbind(1:length(loc[,1]),as.numeric(response))
np = length(loc[,1]) # it must be initialize in Whilelm function
# set lambda

setwd(gamfemfitfolder)
#setwd("/Users/giuliopn/PACSworkspace3/PACSworkspace/GAM_tests/")

# Whilelm.fit = gam.fem.fit( data = observed.data ,desmat=desmat ,fdobj = simfd ,lambda = lambda_W, max.steps=15, 
#                            fam="binomial", mesh.const=T, method.phi=1, psi=NULL, 
#                            scale.param=NULL, GCV.score=FALSE, tune=1.8, weight=F )

Whilelm.fit = gam.fem( data = observed.data ,desmat=desmat ,fdobj = simfd ,log.lambda.seq = lambda, max.steps=15, 
                       fam="binomial", mesh.const=T, method.phi=1, psi=NULL, show=F,
                       scale.param=NULL, tune=1.8 )


# Check betas -------------------------------------------------------------
Whilelm.fit$beta
betas_hat_CPP

# Check f -----------------------------------------------------------------
Whilelm.fit$felsplobj$coefs
output_CPP$fit.FEM$coeff




