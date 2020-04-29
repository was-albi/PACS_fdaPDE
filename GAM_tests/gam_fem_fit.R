gam.fem.fit= function(data ,desmat ,fdobj ,lambda, max.steps=15, mu0=NULL, fam=c("binomial", "probit","cloglog", "exponential", "gamma", "poisson"), mesh.const=T,location=NULL, method.phi=1, psi=NULL, scale.param=NULL, GCV.score=FALSE, tune=1.8, weight=F )
{ 
  # SMOOTH.FEM.FD Compute a solution for a Spatial Spline problem for the exponential family of distributions 
  # FOR A GIVEN SMOOTHING SMOOTHING PARAMETER. The optimization of the smoothing parameter is not included in this function.
  # See the function gam.fem for the function that (roughly) includes the estimation of the smoothing parameter.
  # Extends the methodology of Sangalli et. al. (2013) to the whole exponential family (analogously to extend LM to GLM).
  # The reference is: Wilhelm, M. & Sangalli, L. M., Generalized spatial regression with differential regularization, (2016)
  # Journal of Statistical Computation and Simulation, 86, 2497–2518, 2016.
  # In the following, the equation numbers refer to the ArXiv version of the paper:
  # https://arxiv.org/abs/1511.02688
  # The general model is: E[y_i] = g(x_k^t \beta + L(f)_i, for i = 1,...,n. See the general formulation, appendix D
  # E[.] denotes the expectation, y_i is the observation, g(.) is a given link function,  x_k^t is a vector of known covariates
  # L(.)_i is a given linear functional of f (the evaluation or the integral over a subdomain). 
  # The parameters to be estimated are: \beta (which is a vector of covariates), f is a function.
  # The scale parameter of the exponential family distributions \phi may alos be estimated. See details in Section 4.5. 
  # Arguments:
  # FELSPLOBJ    a FELspline object.
  # LAMBDA       a scalar smoothing parameter.
  # data         a n-by-2 set of noisy observations of the surface values.  
  #              data[,1] indexes the points at which the 
  #              values in data[,2] were observed. 
  # desmat       Matrix of covariates (for the parametric part, denoted by X in the paper).
  # max.steps    Maximum number of steps run in the PIRLS algorithm, set to 15 by default
  # mu0          Initial value of the mean (natural parameter). There exists a default value for each familiy
  # mesh.const   Boolean. True if the vertices of the mesh coincides with the locations of the points. Then if T, the 
  # method.phi   Way of computing the scale parameter phi. If method.phi ==1, then Equation 4.7 is used. If 
  #              method.phi == 2 is selected, then the last equation of  paragraph 3.3, p.61,
  #              of Marra and Wood (2012, Scand. journ. Stat)
  # fam          String. Denotes the distribution of the data, within the exponential family. 
  # scale.param  If necessary and not supplied, the scale parameter \phi is estimated. See method.phi for details.
  # tune         Parameter \gamma in the paper, Equation (4.9). It is usually set to 1, but may be higher. It gives more weight
  #              to the equivalent degrees of freedom in the computation of the value of the GCV.
  #     Output:
  # FELSPLOBJ  ...  A FD object of the FEM type defined by the coefficient
  #                 vector resulting from smoothing
  # LAPLACEFD  ...  A FD object of the FEM type for the value of the 
  #                 Laplace operator if order == 2, or empty if order == 1
  #
  #
  
#define loglike
  log.like = NA
  CONSTPREC = 16
#  check arguments
if(!is.vector(fam))
{stop("fam must the name of a distribution")} else
{
	if(length(fam)>1 || ! all(fam%in% c("binomial", "probit","cloglog", "exponential", "gamma", "poisson")) )
	{stop("invalid parameter fam")}
}
if (!is.fd(fdobj))
{   stop('FDOBJ is not a FD object')
}
if (!is.numeric(lambda))
{  
   stop('LAMBDA is not numeric')
}  else if (length(lambda) != 1)
{   
   stop('LAMBDA is not a scalar')
}

#  check data argument

if (is.null(data))
{data=getdata(fdobj)} else 
{
	if (length(dim(data))>0)
    {
    	if (dim(data)[[2]]!=2)
       {
       	if (dim(data)[[1]]!=2)
        {stop('DATA is not a n-by-2 array')}   else
        {data=t(data)}
       }
     } else {stop('DATA is not a n-by-2 array')}
}
# if(!mesh.const)
# {
	# if(is.null(location))
	# {
		# warning("locations necessary if the mesh is constrained by the location of the observations")
	# } else
	# {
		# if(dim(location)[2]!=2)
		# {stop("location of the observations are not valid points")}
	# }
# }
if(!mesh.const & !is.null(location))
{
	if(dim(location)[2]!=2)
	stop("location of the observations are not valid points")
}
# number of data
n<-dim(data)[1]
#basic function to create a canonical basis vector
basis.vect<-function(j,m)
{
	coef<- numeric(m)
	coef[j]<-1
	return(coef)
}  

#  Construct penalty matrix and 'b' vector for Ax=b.


basisobj = fdobj$basis
numnodes = dim(basisobj$params$nodes)[[1]]
nodeStruct = basisobj$params


#  ---------------------------------------------------------------
# construct mass matrix K0 (inner product of the basis functions)
#  ---------------------------------------------------------------
K0 = mass(nodeStruct)
file_name = "debugging_output/R/R0.txt"
write.table(format(K0,scientific = TRUE), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
#  ---------------------------------------------------------------
# construct stiffness matrix K1 (inner product of the gradients of basis functions)
#  ---------------------------------------------------------------
K1 = stiff1(nodeStruct)
file_name = "debugging_output/R/R1.txt"
write.table(format(K1,digits=CONSTPREC), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
# dimension of the FEM space
m  = basisobj$nbasis
#  ---------------------------------------------------------------
# construct the penalty matrix P with ones on diagonal at data points
#  ---------------------------------------------------------------
if(fam=="gamma" & is.null(scale.param))
{scale.param.flag<-T} else
{scale.param.flag<-F}
if (fam=="poisson")
{
	l<-make.link("log")
	link<-l$linkfun
	inv.link<-l$linkinv
	deriv.link<- function(x){1/x}
	var.link<-function(x){x}
	scale.param=1
	if(is.null(mu0))
	{
		mu0=(data[,2])
		mu0[mu0==0]<-1
	}
	loglike.fun<-function(x,mu, scale.param){dpois(x,lambda=mu, log=T)}
	dev.fun<-function(x,mu,scale.param)
	{{poisson()$dev.resids(y= x, mu=mu, wt=rep(1,length(mu)))}}
}

if (fam=="binomial")
{
	link<-logit
	inv.link<-inv.logit
	deriv.link<- function(x){1/(x*(1-x))}
	var.link<-function(x){x*(1-x)}
	scale.param=1
	if(is.null(mu0))
	{mu0=(data[,2]+0.5)/2}
	loglike.fun<-function(x,mu, scale.param){dbinom(x,size=1,prob=mu, log=T)}
	dev.fun<-function(x,mu,scale.param)
	{
		dev=numeric(length(x))
		dev[x==0]= 2*log(1/(1-mu[x==0]))
		dev[x==1]= 2*log(1/mu[x==1])
		dev
	}
}
if(fam=="exponential")
{	
	link<-function(x){-1/x}
	inv.link<-link
	deriv.link<- function(x){(1/x)^2}
	var.link<-function(x){x^2}
	#since mean(exp(lambda))=1/lambda, as first estimate of lambda, we set 1/y_i
	if(is.null(mu0))
	{mu0=data[,2]}
	loglike.fun<-function(x,lambda=1/mu, scale.param){return(dexp(x, rate=lambda, log=T))}
	dev.fun<-function(x, mu, scale.param){2*(((x-mu)/mu)*log(x/mu) )}

}
if(fam=="gamma")
{	
	link<-function(x){-1/x}
	inv.link<-link
	deriv.link<- function(x){(1/x)^2}
	var.link<-function(x){x^2}
	if(is.null(mu0))
	{mu0=data[,2]}
	loglike.fun<-function(x, mu, scale.param){return(dgamma(x, shape= mu/scale.param,scale=scale.param , log=T))}
	dev.fun<-function(x, mu, scale.param){2*(((x-mu)/mu)*log(x/mu) )}
}
if (fam=="probit")
{
	l<-make.link("probit")
	link<-l$linkfun
	inv.link<-l$linkinv
	deriv.link<- function(x)exp((x^2)/2)
	var.link<-function(x){x*(1-x)}
	scale.param=1
	if(is.null(mu0))
	{mu0=(data[,2]+0.5)/2}
	loglike.fun<-function(x,mu, scale.param){dbinom(x,size=1,prob=mu, log=T)}
	dev.fun<-function(x,mu,scale.param)
	{
		dev=numeric(length(x))
		dev[x==0]= 2*log(1/(1-mu[x==0]))
		dev[x==1]= 2*log(1/mu[x==1])
		dev
	}
}
if (fam=="cloglog")
{
	warning("clolog not ready yet, deriv.link false!!")
	l<-make.link("cloglog")
	link<-l$linkfun
	inv.link<-l$linkinv
	deriv.link<- l$mu.eta
	var.link<-function(x){x*(1-x)}
	scale.param=1
	if(is.null(mu0))
	{mu0=(data[,2]+0.5)/2}
	loglike.fun<-function(x,mu, scale.param){dbinom(x,size=1,prob=mu, log=T)}
	dev.fun<-function(x,mu,scale.param)
	{
		dev=numeric(length(x))
		dev[x==0]= 2*log(1/(1-mu[x==0]))
		dev[x==1]= 2*log(1/mu[x==1])
		dev
	}

}
if(!mesh.const & is.null(psi))
{psi<-matrix(sapply(1:n,function(x) sapply(1:m, function(y) eval.FEM.fd(X=location[x,1], Y=location[x,2], fdobj  = fd(basis.vect(y,m=m), basisobj) ))), nrow=n, ncol=m, byrow=T)}
#number of covariates
q=dim(desmat)[2]
#browser()
#computing the matrices for the variance estimation
if((is.null(scale.param) & (fam %in% c("gamma"))) | GCV.score)
{
	if(mesh.const)
	{
		psi=matrix(0, nr=n, nc=m)
		diag(psi)= 1
	}
	#big design matrix
	X = cbind(desmat, psi)
	#big penalization matrix
	S=matrix(0, nr=m+q, nc=m+q)
	S[(q+1):(m+q), (q+1):(m+q)]= t(K1)%*% solve(K0) %*% K1	
}
#initialisation of the parameters
mu<-mu0
file_name = "debugging_output/R/mu0.txt"
write.table(format(mu0,scientific = TRUE), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)

estimation.fn<-matrix(0, nr=dim(data)[1] , nc=max.steps)
estimation.beta<-matrix(0, nr=q , nc=max.steps)
estimation.mu<-matrix(0, nr=dim(data)[1] , nc=max.steps)
pen.log.like<-numeric(max.steps)
#browser()
for (j in 1:max.steps)
{
  print(j)
	eta=link(mu)
	W<-diag(1/(deriv.link(mu)^2 * var.link(mu)))
	file_name = "debugging_output/R/WeightMatrix.txt"
	write.table(format(W,digits=CONSTPREC), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
	z<-eta+(data[, 2]-mu)*deriv.link(mu)
	file_name = "debugging_output/R/Z.txt"
	write.table(format(z,digits=CONSTPREC), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
	#sqrt.W<-sqrt(W)
	#desmatprod = solve( t(X.tilde)  %*% X.tilde ) %*% t(X.tilde)
	desmatprod= solve( t(desmat) %*% W %*% desmat) %*% t(desmat) %*% W
	file_name = "debugging_output/R/desmatprod.txt"
	write.table(format(desmatprod,scientific = TRUE), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
	
	H = desmat %*% desmatprod
	file_name = "debugging_output/R/H.txt"
	write.table(format(H,digits=CONSTPREC), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
	Q = t(diag(1,length(data[,1]))-H) %*% W %*% (diag(1,length(data[,1]))-H)
	
	#Q = diag(1,length(data[,1]))-H
	
	file_name = "debugging_output/R/Q.txt"
	write.table(format(Q,digits=CONSTPREC), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
	
	#Q = diag(1,length(data[,1]))-H shouldn't be this one? Does it hold becouse Q is a proj matrix? - Alb
	if(mesh.const)
	{
		L =matrix(0,numnodes,numnodes)
		L[data[,1],data[,1]] = Q
		b            = matrix(numeric(numnodes*2),ncol=1)
		b[data[,1],] = Q %*% z
	} else
	{
			L =matrix(0,numnodes,numnodes)
			L=t(psi) %*% Q %*% psi
			b = matrix(numeric(numnodes*2),ncol=1)
			b[1:numnodes,] = t(psi) %*% Q %*% z
	}
	
	file_name = "debugging_output/R/b.txt"
	write.table(format(b,digits=CONSTPREC), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
# #  ---------------------------------------------------------------
# # construct vector b for system Ax=b
# #  ---------------------------------------------------------------

    AS       = rbind(cbind(L, -lambda*K1), cbind(K1, K0))
    bigsol   = solve(AS,b)
    
    file_name = "debugging_output/R/AS.txt"
    write.table(format(AS,digits=CONSTPREC), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
    file_name = "debugging_output/R/bigsol.txt"
    write.table(format(bigsol,digits=CONSTPREC), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
# #  ---------------------------------------------------------------
# # get the solutions of the problem: function f evaluated on the internal nodes (fnhat),
# # function f evaluated on the whole mesh (u), laplacian of f evaluated on the whole mesh (s) and beta
# #  ---------------------------------------------------------------
	u = bigsol[1:numnodes,]
	s = bigsol[(numnodes+1):(2*numnodes),]

	file_name = "debugging_output/R/solution_f.txt"
	write.table(format(u,digits=CONSTPREC), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
	file_name = "debugging_output/R/solution_lapl_f.txt"
	write.table(format(s,digits=CONSTPREC), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
	
	#functional object (f and laplacian of f) defined by the coefficients
    felsplobj  = fd(u, basisobj)
	laplacefd = fd(s, basisobj)
	if(mesh.const)
	{fnhat = as.vector(bigsol[1:np,])} else
	{fnhat = psi %*% u}
	# fdobj$basis$params$nodes
	  tmp = z-fnhat
	  file_name = "debugging_output/R/z_meno_fnhat.txt"
	  write.table(format(tmp,scientific = TRUE), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
	  
    betahat  = desmatprod %*% (z-fnhat)
    file_name = "debugging_output/R/beta_hat.txt"
    write.table(format(betahat,scientific = TRUE), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
    estimation.beta[, j]<-betahat
    estimation.fn[, j]<-fnhat
      mu       = as.vector(inv.link(desmat %*% betahat + fnhat))
      file_name = "debugging_output/R/mu.txt"
      write.table(format(mu,digits=CONSTPREC), file_name, append = FALSE, sep = ",", dec = ".", row.names = FALSE, col.names = FALSE)
    estimation.mu[, j]=mu
	if(scale.param.flag)
	{
		#big hat matrix: applied to pseudo data gives the estimates \hat{beta} and \hat{f},
		#denoted by A and called influence matrix in Wood (2006)
		bigF= solve((t(X) %*% W %*% X+lambda*S)) %*% t(X) %*% W %*% X
		# here, var(X)= phi*var.link(mu) 
		# if method.phi==2, phi is computed as in Marra & Wood (2012) (p. 61 for the derivation of the value)
		# else phi is computed as in Wood IGAM
		if(method.phi==2)
		{
			big.beta<- as.vector(c(betahat,u))
			res <- z- (desmat %*% betahat) - fnhat
			quad1<- sum(diag(W)* (res^2))
			Fbeta<- as.vector(X %*% ((Fdef %*% big.beta)- big.beta))
			quad2<-sum(diag(W)* (Fbeta)^2)
			tr.FF<-sum(apply(Fdef,c(1,2),function(x) x^2 ))
			phi<-(quad1-quad2)/(n-(2*sum(diag(Fdef)))+ tr.FF)
		} else
		{
			phi <- sum(((data[,2]-mu)^2)/var.link(mu))/(n-sum(diag(bigF)))
		}
		#here the dispersion parameter is computed as the mean of all dispersion parameters
		scale.param = mean((var.link(mu)*phi)/mu)
	}
	# penalized log likelihood, i.e the functional (to be minimized) evaluated at the current estimate
	pen.log.like[j]<-sum(loglike.fun(data[,2],mu, scale.param))- lambda *(t(s)%*% K0 %*% s)
	# truncate the results and stop the iteration if convergence is reached
	if(j>1)
	{
		
		if(round(pen.log.like[j], digits=4)==round(pen.log.like[j-1], digits=4))
		{
			estimation.fn<-estimation.fn[,1:(j-1)]
			estimation.mu<-estimation.mu[,1:(j-1)]
			pen.log.like<-pen.log.like[1:(j-1)]
			log.like=sum(loglike.fun(data[,2],mu, scale.param))
			if(GCV.score)
			{
				if(!scale.param.flag)
				{
					#big hat matrix: applied to pseudo da ta gives the estimates \hat{beta} and \hat{f},
					#denoted by A and called influence matrix in Wood (2006)
					#Already computed if the estimation of the scale parameter is included
					bigF = solve((t(X) %*% W %*% X+lambda*S)) %*% t(X) %*% W %*% X
					deg.freed = sum(diag(bigF))
					print("DOFs:\n")
					print(deg.freed)
					#GCV.score = ((1/n)*sum(dev.fun(x=data[,2],mu=mu, scale.param)))-(mean(scale.param*var.link(mu)))+((2/n)*deg.freed*(mean(scale.param*var.link(mu)))*tune)
					GCV.score = (n*sum(dev.fun(x=data[,2],mu=mu, scale.param)))/((n-tune*deg.freed)^2)
				} else
				{
				deg.freed = sum(diag(bigF))
				print("DOFs:\n")
				print(deg.freed)
				#actually URBE score, since scale.param is known
				GCV.score = (n*sum(dev.fun(x=data[,2],mu=mu, scale.param)))/((n-tune*deg.freed)^2)
				}
				
			}
			iteration<-j-1
			break
		}
	}
	if(j==max.steps)
	{
		warning("max iteration reached")
		iteration=max.steps
		#tuning constant  to ensure the GCV convexity w.r.t smoothing parameter
		if(GCV.score)
		{
			if(!scale.param.flag)
			{
				#big hat matrix: applied to pseudo data gives the estimates \hat{beta} and \hat{f},
				#denoted by A and called influence matrix in Wood (2006)
				#Already computed if the estimation of the scale parameter is included
				bigF = solve((t(X) %*% W %*% X+lambda*S)) %*% t(X) %*% W %*% X
				deg.freed = sum(diag(bigF))
				print("DOFs:\n")
				print(deg.freed)
				#actually URBE score, since scale.param is known
				#GCV.score = ((1/n)*sum(dev.fun(x=data[,2],mu=mu, scale.param)))-(mean(scale.param*var.link(mu)))+((2/n)*deg.freed*(mean(scale.param*var.link(mu)))*tune)
				GCV.score = (n*sum(dev.fun(x=data[,2],mu=mu, scale.param)))/((n-tune*deg.freed)^2)
			} else
			{
				deg.freed = sum(diag(bigF))			
				print("DOFs:\n")
				print(deg.freed)
				GCV.score = (n*sum(dev.fun(x=data[,2],mu=mu, scale.param)))/((n-tune*deg.freed)^2)
			}
		}
	}
}
#browser()
# Make FELspline object
#function f defined on the mesh   
felsplobj  = fd(u, basisobj)
laplacefd = fd(s, basisobj)
if(weight==T)
{reslist=list(felsplobj=felsplobj,laplacefd=laplacefd,mu=mu, beta= betahat, estimation.mu=estimation.mu, pen.log.like=pen.log.like,log.like=log.like, estimation.fn=estimation.fn, numb.iterations= iteration, lambda=lambda, scale.param=scale.param, W=W)} else
{reslist=list(felsplobj=felsplobj,laplacefd=laplacefd,mu=mu, beta= betahat, estimation.mu=estimation.mu, pen.log.like=pen.log.like,log.like=log.like, estimation.fn=estimation.fn, numb.iterations= iteration, lambda=lambda, scale.param=scale.param)}
if(scale.param.flag)
{reslist=c(reslist, phi=phi)}
if(!all(GCV.score==FALSE))
{reslist=c(reslist, GCV.score=GCV.score, deg.freed=deg.freed)}

return(reslist)
}
