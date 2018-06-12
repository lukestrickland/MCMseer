#####LBA transform, random, and likelihood functions
#for unit tests.

transform.lbaN23 <- function(par.df)
{
  par.df$b <- par.df$B+par.df$A

  list(A=par.df$A,b=par.df$b,t0=par.df$t0,mean_v=par.df$mean_v,sd_v=par.df$sd_v,
       st0=par.df$st0,N=par.df$N)
}

random.lbaN23 <- function(n,p.df,model)
{
  rlba.norm(n,
            A=p.df$A[1:p.df$N[1]],
            b=p.df$b[1:p.df$N[1]],
            mean_v=p.df$mean_v[1:p.df$N[1]],
            sd_v=p.df$sd_v[1:p.df$N[1]],
            t0=p.df$t0[1],
            st0=p.df$st0[1],
            posdrift = attr(model,"posdrift"))
}


likelihood.lbaN23 <- function(p.vector,data,min.like=1e-10)
  # Returns vector of likelihoods for each RT in data (in same order)
{

  p.list <- p.list.dmc(p.vector,model=attributes(data)$model,n1order=TRUE,
                       cells=attributes(data)$cells,
                       n1.index=attr(data,"n1.index")
  )

  likelihood <- numeric(dim(data)[1])
  is2 <- p.list$N[,1]==2
  likelihood[is2] <- n1PDFfixedt0.norm(
    dt=data$RT[is2]-p.list$t0[is2,1],
    A=p.list$A[is2,1:2],
    b=p.list$b[is2,1:2],
    mean_v=p.list$mean_v[is2,1:2],
    sd_v=p.list$sd_v[is2,1:2],
    posdrift=attr(attr(data,"model"),"posdrift"))
  likelihood[!is2] <- n1PDFfixedt0.norm(
    dt=data$RT[!is2]-p.list$t0[!is2,1],
    A=p.list$A[!is2,1:3],
    b=p.list$b[!is2,1:3],
    mean_v=p.list$mean_v[!is2,1:3],
    sd_v=p.list$sd_v[!is2,1:3],
    posdrift=attr(attr(data,"model"),"posdrift"))

  pmax(likelihood,min.like)
}
