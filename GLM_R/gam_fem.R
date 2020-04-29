gam.fem<-function(data, desmat,fdobj, max.steps=15, mu0=NULL, fam=c("binomial", "probit","cloglog", "exponential", "gamma"), mesh.const=T, method.phi=1, scale.param=NULL, log.lambda.seq=10^(-2:5), show=T, tune=1.8, psi=NULL)
{
	# Uses a two scale grid (initialized on the log scale) for the optimization of the smoothing parameter.
  # So it is quite a rough optimization method 
  # 
  # log.lambda.seq    Sequence of values of the smoothing parameter for which the model is fitted
  # show              Boolean. If true, the value of the GCV criterion is plotted against the smoothing param
  #time=proc.time()
	test.seq<-lapply(log.lambda.seq,function(x){gam.fem.fit(data=data, desmat=desmat, fdobj=fdobj, lambda=x, max.steps=max.steps, fam=fam, mesh.const=mesh.const, psi=psi, scale.param=scale.param, GCV.score=T, tune=tune)})
	#print(proc.time()-time)
	GCV.log.seq<-as.vector(unlist(lapply(test.seq, function(x)x$GCV)))
	if(show)
	{plot(x=log(log.lambda.seq, base=10), y=GCV.log.seq, xlab="log lambda", ylab="GCV score", main="GCV score")}
	ind.log<-which.min(GCV.log.seq)
	if(ind.log==length(log.lambda.seq)| ind.log==1)
	{
		if(ind.log==length(log.lambda.seq))
		{
			step<-log.lambda.seq[ind.log-1]
			lambda.seq<-seq(from=log.lambda.seq[ind.log], to=log.lambda.seq[ind.log-1], by=-step)
			test<-gam.fem.fit(data=data, desmat=desmat, fdobj=fdobj, lambda=lambda.seq[2], max.steps=max.steps, fam=fam, mesh.const=mesh.const, psi=psi, scale.param=scale.param, GCV.score=T, tune=tune)
			if(test$GCV>=GCV.log.seq[ind.log])
			{
				warning("minimum GCV at the limits")
				plot(x=log(log.lambda.seq, base=10), y=GCV.log.seq)
				return(test=test.seq[[ind.log]])
			} else
			{
				test.seq<-c(test.seq[ind.log],list(test),lapply(lambda.seq[3:9],function(x){gam.fem.fit(data=data, desmat=desmat, fdobj=fdobj, lambda=x, max.steps=max.steps, fam=fam, mesh.const=mesh.const, psi=psi, scale.param=scale.param, GCV.score=T, tune=tune)}), test.seq[ind.log-1])
				GCV.seq<-as.vector(unlist(lapply(test.seq, function(x)x$GCV)))
				lambda.seq<-as.vector(unlist(lapply(test.seq, function(x)x$lambda)))
				test<-test.seq[[which.min(GCV.seq)]]
			}
		} else
		{
			step<-log.lambda.seq[ind.log]
			lambda.seq<-seq(from=log.lambda.seq[ind.log], to=log.lambda.seq[ind.log+1], by=step)
			test<-gam.fem.fit(data=data, desmat=desmat, fdobj=fdobj, lambda=lambda.seq[2], max.steps=max.steps, fam=fam, mesh.const=mesh.const, psi=psi, scale.param=scale.param, GCV.score=T, tune=tune)
			if(test$GCV>=GCV.log.seq[ind.log])
			{
				warning("minimum GCV at the limits")
				plot(x=log(log.lambda.seq, base=10), y=GCV.log.seq)
				return(test=test.seq[[ind.log]])
			} else
			{
				test.seq<-c(test.seq[ind.log],list(test),lapply(lambda.seq[3:9],function(x){gam.fem.fit(data=data, desmat=desmat, fdobj=fdobj, lambda=x, max.steps=max.steps, fam=fam, mesh.const=mesh.const, psi=psi, scale.param=scale.param, GCV.score=T, tune=tune)}), test.seq[ind.log+1])
				GCV.seq<-as.vector(unlist(lapply(test.seq, function(x)x$GCV)))
				lambda.seq<-as.vector(unlist(lapply(test.seq, function(x)x$lambda)))
				test<-test.seq[[which.min(GCV.seq)]]
			}
		}
	} else
	{
		step1<-log.lambda.seq[ind.log-1]
		step2<-log.lambda.seq[ind.log]

		test1<-gam.fem.fit(data=data, desmat=desmat, fdobj=fdobj, lambda=log.lambda.seq[ind.log]-step1, max.steps=max.steps, fam=fam, mesh.const=mesh.const, psi=psi, scale.param=scale.param, GCV.score=T, tune=tune)
		test2<-gam.fem.fit(data=data, desmat=desmat, fdobj=fdobj, lambda=log.lambda.seq[ind.log]+step2, max.steps=max.steps, fam=fam, mesh.const=mesh.const, psi=psi, scale.param=scale.param, GCV.score=T, tune=tune)
		if(test1$GCV>GCV.log.seq[ind.log] & test2$GCV>GCV.log.seq[ind.log])
		{test=test.seq[[ind.log]]}
		if(test1$GCV<GCV.log.seq[ind.log] & GCV.log.seq[ind.log] < test2$GCV)
		{
			step=step1
			lambda.seq<-seq(from=log.lambda.seq[ind.log], to=log.lambda.seq[ind.log-1], by=-step)
			test.seq<-c(test.seq[ind.log],list(test1),lapply(lambda.seq[3:9],function(x){gam.fem.fit(data=data, desmat=desmat, fdobj=fdobj, lambda=x, max.steps=max.steps, fam=fam, mesh.const=mesh.const, psi=psi, scale.param=scale.param, GCV.score=T, tune=tune)}), test.seq[ind.log-1])
			GCV.seq<-as.vector(unlist(lapply(test.seq, function(x)x$GCV)))
			lambda.seq<-as.vector(unlist(lapply(test.seq, function(x)x$lambda)))
			test<-test.seq[[which.min(GCV.seq)]]
		}
		if(test1$GCV>GCV.log.seq[ind.log] & GCV.log.seq[ind.log] > test2$GCV)
		{
			step=step2
			lambda.seq<-seq(from=log.lambda.seq[ind.log], to=log.lambda.seq[ind.log+1], by=step)
			test.seq<-c(test.seq[ind.log],list(test2),lapply(lambda.seq[3:9],function(x){gam.fem.fit(data=data, desmat=desmat, fdobj=fdobj, lambda=x, max.steps=max.steps, fam=fam, mesh.const=mesh.const, psi=psi, scale.param=scale.param, GCV.score=T, tune=tune)}), test.seq[ind.log+1])
			GCV.seq<-as.vector(unlist(lapply(test.seq, function(x)x$GCV)))
			lambda.seq<-as.vector(unlist(lapply(test.seq, function(x)x$lambda)))
			test<-test.seq[[which.min(GCV.seq)]]

		}

	}
	return(test)
}
