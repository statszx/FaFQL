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
covZ_mat=covZY_mat=NULL
for (i in 1:N){
resi_Y[,i]=respon_Y[,i]-init_F%*%matrix(init_lambda[i,])
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


# identify the factor number 

fact_BIC=NULL
for (r in 1:fact_numax){
fact_res=GetFaFQLEst=function(respon_Y,fact_num,cov_Z,pena_mat,tau,h)
fact_theta=fact_res$iter_theta
fact_F=fact_res$iter_F
fact_lambda=fact_res$iter_lambda

fact_resYthe=NULL
for (i in 1:N){
fact_resYthe[,i]=respon_Y[,i]-cov_Z[,,i]%*%fact_theta
}

fact_Vval=eigen(fact_resYthe %*% t(fact_resYthe)/(NT))$val
fact_eignum=r+1
fact_eivnum=T-r
fact_Vmat=matrix(diag(fact_Vval[fact_eignum:T]),ncol=fact_eivnum)
fact_Fmat=eigen(fact_resYthe %*% t(fact_resYthe))$vec
fact_F=matrix(fact_Fmat[,fact_eignum:T])


fact_V=sum(fact_Vval[fact_eignum:T]/(NT))
rho=(N+T)*(dim(W)[2]+dim(X)[2])*ln((NT)/(N+T))/(NT)
fact_BIC[r]=ln(fact_V)+rho*r
}

fact_numopt=which.min(fact_BIC)






