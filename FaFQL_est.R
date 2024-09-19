library(Matrix)
library(MASS)
library(mgcv)
library(psych)
library(pls)
library(plsdof)
library(quantreg)
library(fda)
library(splines2)


# calculate the initial FaFQL estimates

init_res=rq(respon_Y~cov_Z,tau)
init_theta=init_res$coefficients
init_alpha=init_res$coefficients[1]
init_gamma=init_res$coefficients[-1]

init_Vval=eigen(respon_Y %*% t(respon_Y))$val
fact_num=dim(F)[2]
init_Vmat=matrix(diag(init_Vval[1:fact_num]),ncol=fact_num)
init_Fmat=eigen(respon_Y %*% t(respon_Y)/(NT))$vec
init_F=matrix(init_Fmat[,1:fact_num])
init_lambda=t(t(init_F)%*%respon_Y/N)

resi_Y=cov_Z=NULL
covZ_mat=covZY_mat=NULL
for (i in 1:N){
resi_Y[,i]=respon_Y[,i]-init_F%*%matrix(init_lambda[i,])
cov_Z[,,i]=cbind(cov_W[,,i],cov_B[,,i])
covZ_mat=covZ_mat+t(cov_Z[,,i])%*%cov_Z[,,i]
covZY_mat=covZY_mat+t(cov_Z[,,i])%*%resi_Y[,i]
}

init_theta=GetFaFQLEst(respon_Y,fact_num,cov_Z,pena_mat,tau,h)


# calculate the FaFQL estimates from the iteration algorithm

iter_theta=init_theta
iter_num=1
repeat{

iterep_theta=iter_theta
iter_resY=NULL
for(i in 1:N){
iter_resYF[,i]=respon_Y[,i]-cov_Z[,,i]%*%iter_theta
}

iter_Vval=eigen(iter_YF %*% t(resi_YF))$val
iter_Vmat=matrix(diag(iter_Vval[1:fact_num]),ncol=fact_num)
iter_Fmat=eigen(iter_resYF %*% t(iter_resYF)/(NT))$vec
iter_F=matrix(iter_Fmat[,1:fact_num])
iter_lambda=t(t(iter_YF)%*%iter_resYF/N)

iter_resYthe=NULL
itcov_Z=cov_Z
itcovZ_mat=covZ_mat
itcovZY_mat=NULL
for (i in 1:N){
iter_resYthe[,i]=respon_Y[,i]-iter_F%*%matrix(iter_lambda[i,])
itcovZY_mat=itcovZY_mat+t(itcov_Z[,,i])%*%iter_resYthe[,i]
}

iter_theta=ginv(itcovZ_mat+pena_mat)%*%itcovZY_mat

iter_num=iter_num+1
itnor_theta=norm(iter_theta-iterep_theta,2)
if(iter_num>500 | itnor_theta<0.001 )  break

}


# define B-spline basis functions

domain=range(s_grid)
M=30
d=6
norder=d+1
knots=seq(domain[1],domain[2], length.out=M+1)
nknots=length(knots)
nbasis=length(knots)+norder-2
beta_basis =create.bspline.basis(knots,nbasis,norder)
basis_rng=getbasisrange(beta.basis)
breaks=c(basis_rng[1],beta_basis$params,basis_rng[2])
basis_M=dim(X)[2]
basis_s=seq(basis_rng[1],basis_rng[2],length.out=basis_M)
beta_estbasis=eval.basis(basis_s,beta_basis)
L=beta_basis$nbasis
basis_pena=eval.penalty(beta_basis,int2Lfd(2))

GetIntMat=function(cov_indX,beta.basis){
basis_rng=getbasisrange(beta.basis)
breaks=c(basis_rng[1],beta_basis$params,basis_rng[2])
basis_s=seq(basis_rng[1],basis_rng[2],length.out=basis_M)
beta_estbasis=eval.basis(basis_s,beta_basis)
bet_indxmat=NULL
for(i in i:basis_M){
bet_indxmat=bet_indxmat+krnocker(cov_indX[,i],t(beta_estbasis[,i]))
}
betx_mat=bet_indmat/basis_M
return(betx_mat)
}

# create the penalty matrix after B-spline expansion 

basis_paramat=matrix(diag(basis_para),ncol=length(basis_para))
basis_penamat=matrix(basis_pena,ncol=length(basis_pena))
beta_pena=krnocker(basis_paramat,basis_penamat)
alpha_pena=matrix(0,ncol=dim(cov_W)[1],nrow=dim(cov_W)[1])
pena_mat=rbind(cbind(alpha_pena,matrix(0,nrow=nrow(alpha_pena),ncol=ncol(beta_pena))),
cbind(matrix(0,nrow=nrow(beta_pena),ncol=ncol(alpha_pena)),beta_pena))


# objective function 

GetFaFQLEst=function(respon_Y,fact_num,cov_Z,pena_mat,tau,h){

init_res=rq(respon_Y~cov_Z,tau)
init_theta=init_res$coefficients
init_alpha=init_res$coefficients[1]
init_gamma=init_res$coefficients[-1]

resinit_Y=respon_Y-t(cov_Z)%*%init_theta

init_Vval=eigen(resinit_Y %*% t(resinit_Y))$val
init_Vmat=matrix(diag(init_Vval[1:fact_num]),ncol=fact_num)
init_Fmat=eigen(resinit_Y %*% t(resinit_Y))$vec
init_F=matrix(init_Fmat[,1:fact_num])
init_lambda=t(t(init_F)%*%resinit_Y/N)

resi_Y=NULL
smorho=NULL
covZ_mat=covZY_mat=NULL
for (i in 1:N){
resi_Y[,i]=resinit_Y[,i]-init_F%*%matrix(init_lambda[i,])
covZ_mat=covZ_mat+t(cov_Z[,,i])%*%cov_Z[,,i]
covZY_mat=covZY_mat+t(cov_Z[,,i])%*%resi_Y[,i]
smorho=smorho+GetSmoRho(resi_Y[,i],tau,h)
}

iter_theta=init_theta
iter_num=1

repeat{

iterep_theta=iter_theta
iter_resY=NULL
for(i in 1:N){
iter_resYF[,i]=resinit_Y[,i]-cov_Z[,,i]%*%iter_theta
}

iter_Vval=eigen(iter_YF %*% t(resi_YF))$val
iter_Vmat=matrix(diag(iter_Vval[1:fact_num]),ncol=fact_num)
iter_Fmat=eigen(iter_resYF %*% t(iter_resYF)/(NT))$vec
iter_F=matrix(iter_Fmat[,1:fact_num])
iter_lambda=t(t(iter_YF)%*%iter_resYF/N)

iter_resYthe=NULL
itcov_Z=cov_Z
itcovZ_mat=covZ_mat
itcovZY_mat=NULL
smorho=NULL

for (i in 1:N){
iter_resYthe[,i]=respon_Y[,i]-iter_F%*%matrix(iter_lambda[i,])
itcovZY_mat=itcovZY_mat+t(itcov_Z[,,i])%*%iter_resYthe[,i]
smorho=smorho+GetSmoRho(iter_resYthe[,i],tau,h)
}

fafql_res=smorho+t(iter_theta)%*%pena_mat%*%iter_theta*pena_para)

fafql_ge=fafql_e=NULL
for (i in 1:n){
fafql_ge=fafql_ge+GetGauKer(respon_Y[i]-sum(cov_Z[i,]*iter_theta),h)*
as.matrix(c(1,cov_Z[i,]))%*%c(1,cov_Z[i,])
fafql_e=fafql_e-(GauKer(respon_Y[i]-sum(cov_Z[i,]*iter_theta),h)+tau-1)*
c(1,cov_Z[i,])}
fafql_e=as.matrix(fafql_e)
fafql_matrix=bdiag(matrix(0,nrow=nrow(cov_Z),ncol=ncol(cov_Z)+1),
pena_para*pena_mat)+fafql_ge
fafql_mat=svd(fafql_matrix)
fafql_matv=(fafql_mat$u)%*%diag((fafql_mat$d)^-1)%*%t(fafql_mat$v)
fafql_matve=fafql_matv%*%fafql_e
iter_theta=fafql_matve

iter_num=iter_num+1
itnor_theta=norm(iter_theta-iterep_theta,2)
if(iter_num>500 | itnor_theta<0.001 )  break

}
 out=list(iter_theta,iter_F,iter_lambda)
 return(out)
}

# quantile loss function

GetRho=function(u,tau){
  u*(tau-(u<=0))
}

GetLg=function(u){
  (2/pi)^0.5*exp(-u^2/2)+ u*(1-2*pnorm(-u))
}

# Gaussian kernel with bandwidth h

GetGauKer=function(x,h){
  dnorm(x/h)
}

GauKer=function(x,h){
   pnorm(x/h)
}

# smoothed quantile loss function

GetSmoRho=function(u,tau,h){
  h/2*GetLg(u/h)+(tau-0.5)*u
}

