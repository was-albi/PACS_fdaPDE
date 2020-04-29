# Example

#time<-proc.time()
setwd("/home/alb/Scrivania/PACS/Codes_GAM_FEM")  
source("2013_SSR_AllFunctions.R")
library(fda)
library(rgl)
library(mgcv)
library(purrr)
library(boot)
Cconstredges=read.table("multiple_meshes//boundaryedges_TimWood.txt", header=F)
colnames(Cconstredges)=c("x","y")
# source("/Users/Matthieu/Documents/EPFL/PDM/codes/codes_Matth/smooth_GLM_GCV_test.R")
# source("gam_fem.R")
# source('gam_fem_fit.R')
source("C-shaped_domain_fun.R")

# mu 1:10 -> param -1:-1/10

# muvec -1:-0.35
# beta1 * cov1  -0.8 : -0.4  -2/5 * beta(alpha,beta)
# beta2 * cov2  0.3 : 0.6  3/10 (beta(alpha,beta)+0.3)
np = 200

link<-function(x){logit(x)}
inv.link<-function(x){1/(1+exp(-x))}



######################################################################
######################################################################
######################   DATA GENERATION   ###########################
######################################################################
######################################################################
fsb <- list(fs.boundary())
#fsb<-list(fs.boundary(r0=0.3, r=0.5, l=3, n.theta=10))
nx<-250
ny<-100 
xvec <- seq(-1,4,length=nx)
yvec<-seq(-1,1,length=ny)
xx <- rep(xvec,ny)
yy<-rep(yvec,rep(nx,ny))
a = 8;
b=10 # 10
tru <- -(1/a)*(matrix(fs.test(xx,yy),nx,ny) +b) ## truth

image(xvec,yvec,tru,col=heat.colors(100),xlab="x",ylab="y")
lines(fsb$x,fsb$y,lwd=3)
contour(xvec,yvec,tru,levels=seq(-5,5,by=.25),add=TRUE)

scale.param=.1 # var = mean * sigmaquad

knots <- data.frame(x=rep(seq(-.5,3,by=.5),4),y=rep(c(-.6,-.3,.3,.6),rep(8,4)))

nrep=100

covariates <- array(0, dim=c(np, 2, nrep))
mu.ex <- matrix(0, nc=np, nr=nrep)
region <- matrix(0, nc=np, nr=nrep)
observed.response <- matrix(0,nr=nrep ,ncol=np)


fit.fem <- list()
space.est.fem <- array(0, dim=c(nx,ny,nrep))
beta.est.fem <- matrix(0,nr=nrep ,ncol=2)
fnhat.est.fem <- matrix(0, nc=np, nr=nrep)
mu.est.fem <- matrix(0, nc=np, nr=nrep)
smooth.fem <- numeric(length=nrep)
scale.fem <- numeric(length=nrep)
count.check <- numeric(length=nrep)

space.est.soap<- array(0, dim=c(nx,ny,nrep))
beta.est.soap<-matrix(0,nr=nrep ,ncol=2)
mu.est.soap<-matrix(0, nc=np, nr=nrep)
fitted.soap<-matrix(0, nc=np, nr=nrep)
smooth.soap<-matrix(0, nc=2, nr=nrep)
soap.converged<-numeric(length=nrep)
soap.coefficients <- matrix(0,nr=nrep,nc=73 ) # 73 is the number of coefficients
Matpred.soap <- array(0, dim=c(nx*ny, 71, nrep))

space.est.tps<- array(0, dim=c(nx,ny,nrep))
beta.est.tps<-matrix(0,nr=nrep ,ncol=2)
mu.est.tps<-matrix(0, nc=np, nr=nrep)
fitted.tps<-matrix(0, nc=np, nr=nrep)
smooth.tps<-matrix(0, nc=2, nr=nrep)
tps.converged<-numeric(length=nrep)
tps.coefficients <- matrix(0,nr=nrep,nc=42 ) # 73 is the number of coefficients
Matpred.tps <- array(0, dim=c(nx*ny, 40, nrep))
set.seed(1)
seeds <- 10 * sample(2*nrep, size = nrep, replace = F)
seeds.check <- rep(NA, nrep)







sim=1


loc.name=paste("multiple_meshes/","location_",sim,".txt", sep="" )
tri.name=paste("multiple_meshes/","tri_",sim,".txt", sep="" )
pmesh        = read.table(loc.name, header=T)
Tri = as.matrix(read.table(tri.name, header=F))
names(pmesh)=c("x","y")
x=pmesh[,1]
y=pmesh[,2]
p            = pmesh[,1:2]
order = 1
p2col        = rbind(p,Cconstredges)
e = NULL
#  set up the FEM basis object and plot it
basisobj = create.FEM.basis(p2col, e, Tri, order)
#  set up a dummy FEM functional data object
simfd = fd(numeric(basisobj$nbasis),basisobj)


muvec = -(1/a)*( fs.test(x,y)+b) # value of f in the observations
desmat=matrix(0,nrow=np,ncol=2)
mu=numeric(length=np)
count <- 0
# create some suited covariates
while(any(mu<=0))
{
  count <- count+1
  # The seed is set in such a way that if count <10, which is very likely, then
  # the seeds are different for all the simulations (which is highly desirable!)
  set.seed(seeds[sim] + count)
  desmat[,1]=rbeta(np,shape1=1.5,shape2=2)  # sampling covariates from beta distr.
  desmat[,2]=rbeta(np,shape1=3,shape2=2)+1  # sampling covariates from beta distr.
  beta1=-2/5
  beta2=3/10
  covariates[,,sim]<-desmat
  param=muvec+beta1*desmat[,1]+beta2*desmat[,2]
  mu<-inv.link(param)
  if(all(mu>=0)) seeds.check[sim] <- seeds[sim] + count
}
count.check[sim] <- count
rm(count)
set.seed(seeds[sim])
response <- rbernoulli(np,p = mu)
#response<-rgamma(np, shape=mu/scale.param, scale=scale.param)
observed.response[sim,] <- response
location = as.data.frame(cbind(x,y), names=c("x","y"))
region[sim,]= apply(location, 1,zone.att)
mu.ex[sim,]<-mu
observed.data = cbind(1:np,response)

####### FEM estimation ######

#test<-gam.fem(data=observed.data, desmat=desmat, fdobj=simfd, max.steps=15, fam="binomial", mesh.const=T, scale.param=NULL, log.lambda.seq=10^(0:6))
test <- gam.fem.fit(data=observed.data, desmat=desmat, fdobj=simfd, lambda = 10^-5, max.steps=15, fam="binomial", mesh.const=T, scale.param=NULL)
fit.fem[[sim]] <- test
space.est.fem[,,sim]<-matrix(eval.FEM.fd(X=xx, Y=yy,test$felsplobj),nx,ny)
beta.est.fem[sim,]<-test$beta
fnhat.est.fem[sim,]<-as.vector(eval.FEM.fd(X=location[,1], Y=location[,2],test$felsplobj))
mu.est.fem[sim,]<-inv.link(as.vector(fnhat.est.fem[sim,]+desmat %*% beta.est.fem[sim,]) )
smooth.fem[sim]<-test$lambda
scale.fem[sim] <-test$scale 

library(fdaPDE)

mesh=create.mesh.2D(nodes=knots)
plot(mesh)
plot(knots)
FEMbasis=fdaPDE::create.FEM.basis(mesh)

test <- gam.fem.fit(observations = observed.data, FEMbasis =FEMbasis, covariates = desmat, lambda = 10^-5, max.steps=15, fam="binomial", mesh.const=T, scale.param=NULL)



















