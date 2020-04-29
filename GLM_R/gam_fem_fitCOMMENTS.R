#################################################################
# ATTENTION:													#
#	"gam_fem_fitCOMMENTS.R is the file with the main function. 	#
#	 It will be used only for comments, 						#
#	 and to understand the implementation of PIRLS algorithms"	#
################################################################# 

# The file are also used to understand the different notation between
# the Willhem paper and PseudoReport.

#Changes:
#	1) Rename the matrices K0,K1 in the right way(R0,R1)
#	2) Rename the matrix W as A 
#	3) Rename z (the pseudo-data) as z_tilde
# 	4) Rename desmat (covariates) as W 
#	5) Make sections for the all function
#	6) Make the basics steps of PIRLS Algorithm
#	7) data[i,2] = z_i observed values (data[,1] are used as indecs )

gam.fem.fit= function(data ,W ,fdobj ,
	lambda, max.steps=15, mu0=NULL, 
	fam=c("binomial", "probit","cloglog", "exponential", "gamma", "poisson"), 
	mesh.const=T,location=NULL, method.phi=1, 
	psi=NULL, scale.param=NULL, GCV.score=FALSE, 
	tune=1.8, weight=F )
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
  #              values in data[,2] were observed. ( in the paper we called z_i vector the data[,2] )
  # W       Matrix of covariates (for the parametric part, denoted by X in the paper, W in the PseudoReport
  # max.steps    Maximum number of steps run in the PIRLS algorithm, set to 15 by default
  # mu0          Initial value of the mean (natural parameter). There exists a default value for each familiy
  # mesh.const   Boolean. True if the vertices of the mesh coincides with the locations of the points. Then if T, the 
  # method.phi   Way of computing the scale parameter phi. If method.phi ==1, then Equation 4.7 is used. If 
  #              method.phi == 2 is selected, then the last equation of paragraph 3.3, p.61,
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
  


### 1) CHECK ARGUMENTS
{
	if(!is.vector(fam))
	{stop("fam must the name of a distribution")} else
	{
		if(length(fam)>1 | ! all(fam%in% c("binomial", "probit","cloglog", "exponential", "gamma", "poisson")) )
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
}
### 2) Construct R0, R1, P, psi, X matrices
### 	and initialize the variabiles: link, inv.link, 
###										deriv.link, var.link, 
###										loglike.fun, and dev.fun functions
###									 	and mu0 based on family_type
{
	#  Construct penalty matrix and 'b' vector for Ax=b.
	basisobj = fdobj$basis
	numnodes = dim(basisobj$params$nodes)[[1]]
	nodeStruct = basisobj$params


	#  ---------------------------------------------------------------
	# construct mass matrix R0 (inner product of the basis functions) (in the paper R0) 
	#  ---------------------------------------------------------------
	R0 = mass(nodeStruct)
	#  ---------------------------------------------------------------
	# construct stiffness matrix R1 (inner product of the gradients of basis functions) (in the paper R1) 
	#  ---------------------------------------------------------------
	R1 = stiff1(nodeStruct)
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
		#since mean(exp(lambda))=1/lambda, as first estimate of lambda, we set 1/z_i
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

	# From the previous "if", the variabiles: link, inv.link, deriv.link, var.link, loglike.fun, and dev.fun are computed

	if(!mesh.const & is.null(psi))
	{psi<-matrix(sapply(1:n,function(x) sapply(1:m, function(y) eval.FEM.fd(X=location[x,1], Y=location[x,2], fdobj  = fd(basis.vect(y,m=m), basisobj) ))), nrow=n, ncol=m, byrow=T)}
	#number of covariates
	q=dim(W)[2]
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
		X = cbind(W, psi)
		#big penalization matrix
		S=matrix(0, nr=m+q, nc=m+q)
		S[(q+1):(m+q), (q+1):(m+q)]= t(R1)%*% solve(R0) %*% R1
	} 
}

### 3) Initialisation of the parameters: fn, beta, mu
{
	mu<-mu0
	estimation.fn<-matrix(0, nr=dim(data)[1] , nc=max.steps)
	estimation.beta<-matrix(0, nr=q , nc=max.steps)
	estimation.mu<-matrix(0, nr=dim(data)[1] , nc=max.steps)
	pen.log.like<-numeric(max.steps)
}
#browser()
### 4) Start PIRLS Algo
{

	for (j in 1:max.steps)
	{
		
	### 1) Initialize eta, A, z_tilde, H, Q, L= psi %*% Q %*% psi
	{	
		# eta = g(mu)
		eta=link(mu)
		# Compute A at j step
		A<-diag(1/(deriv.link(mu)^2 * var.link(mu)))#<--- DUBBIO 
		
		# Compute z_tilde (the pseudo-data) at j step
		z_tilde<-eta+(data[, 2]-mu)*deriv.link(mu)
		#sqrt.A<-sqrt(A)
		# NB: solve(...) = (...)^-1 
		#Wprod = solve( t(X.tilde)  %*% X.tilde ) %*% t(X.tilde)
		#Wprod = ( W^t %*% A %*% W)^-1 W^t %*% A  
		Wprod= solve( t(W) %*% A %*% W) %*% t(W) %*% A
		H = W %*% Wprod
		Q = t(diag(1,length(data[,1]))-H) %*% A %*% (diag(1,length(data[,1]))-H)
	}	
	### 2) Construct the system of Prop 4.2
	{
		if(mesh.const) 
		{
			L =matrix(0,numnodes,numnodes)
			L[data[,1],data[,1]] = Q
			b            = matrix(numeric(numnodes*2),ncol=1)
			b[data[,1],] = Q %*% z_tilde
		} else
		{
				L =matrix(0,numnodes,numnodes)
				L=t(psi) %*% Q %*% psi
				b = matrix(numeric(numnodes*2),ncol=1)
				b[1:numnodes,] = t(psi) %*% Q %*% z_tilde
		}
		# #  ---------------------------------------------------------------
		# # construct vector b for system Ax=b ---> Prop 4.2 see Wilhelm paper
		# #  --------------------------------------------------------------
	    AS       = rbind(cbind(L, -lambda*R1), cbind(R1, R0))
	    bigsol   = solve(AS,b)
		# #  ---------------------------------------------------------------
		# # get the solutions of the problem: function f evaluated on the internal nodes (fnhat),
		# # function f evaluated on the whole mesh (u), laplacian of f evaluated on the whole mesh (s) and beta
		# #  ---------------------------------------------------------------
		u = bigsol[1:numnodes,]
		s = bigsol[(numnodes+1):(2*numnodes),]
		#functional object (f and laplacian of f) defined by the coefficients
	    felsplobj  = fd(u, basisobj)
		laplacefd = fd(s, basisobj)
		if(mesh.const)
		{fnhat = as.vector(bigsol[1:np,])} else
		{fnhat = psi %*% u}
		# fdobj$basis$params$nodes
	    betahat  = Wprod %*% (z_tilde-fnhat)
	    estimation.beta[, j]<-betahat
	    estimation.fn[, j]<-fnhat
	    mu       = as.vector(inv.link(W %*% betahat + fnhat))
	    estimation.mu[, j]=mu
	}
	### 3) COMPUTE the scale parameter \phi (iff method.phi == 1 )
	{	
		if(scale.param.flag)
		{
			#big hat matrix: applied to pseudo data gives the estimates \hat{beta} and \hat{f},
			#denoted by A and called influence matrix in Wood (2006)
			bigF= solve((t(X) %*% A %*% X+lambda*S)) %*% t(X) %*% A %*% X
			# here, var(X)= phi*var.link(mu) 
			# if method.phi==2, phi is computed as in Marra & Wood (2012) (p. 61 for the derivation of the value)
			# else phi is computed as in Wood IGAM
			if(method.phi==2)
			{
				big.beta<- as.vector(c(betahat,u))
				res <- z_tilde- (W %*% betahat) - fnhat
				quad1<- sum(diag(A)* (res^2))
				Fbeta<- as.vector(X %*% ((Fdef %*% big.beta)- big.beta))
				quad2<-sum(diag(A)* (Fbeta)^2)
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
		pen.log.like[j]<-sum(loglike.fun(data[,2],mu, scale.param))- lambda *(t(s)%*% R0 %*% s)
	}
	### 4) Truncate the results and stop the iteration if convergence is reached
	{	
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
						#big hat matrix: applied to pseudo data gives the estimates \hat{beta} and \hat{f},
						#denoted by A and called influence matrix in Wood (2006)
						#Already computed if the estimation of the scale parameter is included
						bigF = solve((t(X) %*% A %*% X+lambda*S)) %*% t(X) %*% A %*% X
						deg.freed = sum(diag(bigF))
						#GCV.score = ((1/n)*sum(dev.fun(x=data[,2],mu=mu, scale.param)))-(mean(scale.param*var.link(mu)))+((2/n)*deg.freed*(mean(scale.param*var.link(mu)))*tune)
						GCV.score = (n*sum(dev.fun(x=data[,2],mu=mu, scale.param)))/((n-tune*deg.freed)^2)
					} else
					{
					deg.freed = sum(diag(bigF))
					#actually URBE score, since scale.param is known
					GCV.score = (n*sum(dev.fun(x=data[,2],mu=mu, scale.param)))/((n-tune*deg.freed)^2)
					}
					
				}
				iteration<-j-1
				break
			}
		}
	}	

	### 5) Max Step Reach
	{
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
					bigF = solve((t(X) %*% A %*% X+lambda*S)) %*% t(X) %*% A %*% X
					deg.freed = sum(diag(bigF))
					#actually URBE score, since scale.param is known
					#GCV.score = ((1/n)*sum(dev.fun(x=data[,2],mu=mu, scale.param)))-(mean(scale.param*var.link(mu)))+((2/n)*deg.freed*(mean(scale.param*var.link(mu)))*tune)
					GCV.score = (n*sum(dev.fun(x=data[,2],mu=mu, scale.param)))/((n-tune*deg.freed)^2)
				} else
				{
					deg.freed = sum(diag(bigF))				
					GCV.score = (n*sum(dev.fun(x=data[,2],mu=mu, scale.param)))/((n-tune*deg.freed)^2)
				}
			}
		}
	}
	
	}
}
#browser()

### 5) Make FELspline object
{
	#function f defined on the mesh   
	felsplobj  = fd(u, basisobj)
	laplacefd = fd(s, basisobj)
	if(weight==T)
	{reslist=list(felsplobj=felsplobj,laplacefd=laplacefd,mu=mu, beta= betahat, estimation.mu=estimation.mu, pen.log.like=pen.log.like,log.like=log.like, estimation.fn=estimation.fn, numb.iterations= iteration, lambda=lambda, scale.param=scale.param, A=A)} else
	{reslist=list(felsplobj=felsplobj,laplacefd=laplacefd,mu=mu, beta= betahat, estimation.mu=estimation.mu, pen.log.like=pen.log.like,log.like=log.like, estimation.fn=estimation.fn, numb.iterations= iteration, lambda=lambda, scale.param=scale.param)}
	if(scale.param.flag)
	{reslist=c(reslist, phi=phi)}
	if(!all(GCV.score==FALSE))
	{reslist=c(reslist, GCV.score=GCV.score, deg.freed=deg.freed)}

	return(reslist)
}

}
