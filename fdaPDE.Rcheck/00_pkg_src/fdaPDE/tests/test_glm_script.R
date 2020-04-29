#################################
### TEST GLM PDE using fdaPDE ###
#################################
#---- Load Files ----
library(fdaPDE)
#### square 2D (basic case) ####

rm(list=ls())
graphics.off()

data(square2Ddata)

mesh=create.mesh.2D(nodes=nodes)
# x11()
plot(mesh)
# axis(1)
# axis(2)

nnodes=dim(mesh$nodes)[1]

FEMbasis=create.FEM.basis(mesh)

# Test function

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

image(FEM(sol_exact, FEMbasis))

# Set smoothing parameter

lambda= seq(10^-6,10^-3,by=5*10^-6)

GCVFLAG=T
GCVMETHODFLAG='Exact'

data = sol_exact + rnorm(nnodes, mean=0, sd=0.05*(ran[2]-ran[1]))

output_CPP<-gam.fem.fit(observations = data, FEMbasis = FEMbasis, lambda=lambda, fam = fam)

#image(FEM(output_CPP$fit.FEM$coeff[,which.min(output_CPP$GCV)],FEMbasis))

