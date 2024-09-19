# Stock Price Data Analysis

respon_Y=lg(pricl_dly/priop_dly)
cov_W=priop_dly
cov_X=prisec_dly
cov_Z=data.frame(cov_W,GetIntMat(cov_X,beta.basis))

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

fact_BIC=NULL
for(r in 1:fact_numax){
fat_num=r 
fact_res=GetFaFQLEst(respon_Y,fact_num,cov_Z,pena_mat,tau,h)
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

para_thi=para_SSE=NULL

for(j in 1:length(val_grid)){
para_thi[j,]=rep(val_grid[j],length.out=dim(X)[2])
basis_paramat=matrix(diag(para_thi[j,]),ncol=length(dim(X)[2]))
basis_penamat=matrix(basis_pena,ncol=length(basis_pena))
beta_pena=krnocker(basis_paramat,basis_penamat)
alpha_pena=matrix(0,ncol=dim(cov_W)[1],nrow=dim(cov_W)[1])
pena_mat=rbind(cbind(alpha_pena,matrix(0,nrow=nrow(alpha_pena),ncol=ncol(beta_pena))),
cbind(matrix(0,nrow=nrow(beta_pena),ncol=ncol(alpha_pena)),beta_pena))

para_res=GetFaFQLEst(respon_Y,fact_num,cov_Z,pena_mat)
para_theta=para_res$iter_theta
para_F=para_res$iter_F
para_lambda=para_res$iter_lambda

para_resYthe=para_ss=NULL
for (i in 1:N){
para_resYthe[,i]=respon_Y[,i]-cov_Z[,,i]%*%para_theta
para_ss=para_ss+sum(para_resYthe[,i]^2)
}

para_SSE[j]=para_ss
para_Fpro=martix(diag(rep(1,T)))-para_F%*%t(para_F)/T
para_Z=NULL
for(i in 1:N){
para_Z=rbind(para_Z,cov_Z[,,i])
}
para_Zpro=para_Fpro%*%para_Z
para_sthi=para_Z%*%ginv(t(para_Zpro)%*%para_Zpro+pena_mat)%*%t(para_Zpro)
para_s=martix(diag(rep(1,N*T)))-para_sthi
para_gcv[j]=para_SSE[j]/sum(diag(para_s))
}

basis_para=val_grid[which.min(para_gcv)]
basis_paramat=matrix(diag(basis_para),ncol=length(basis_para))
basis_penamat=matrix(basis_pena,ncol=length(basis_pena))
beta_pena=krnocker(basis_paramat,basis_penamat)
alpha_pena=matrix(0,ncol=dim(cov_W)[1],nrow=dim(cov_W)[1])
pena_mat=rbind(cbind(alpha_pena,matrix(0,nrow=nrow(alpha_pena),ncol=ncol(beta_pena))),
cbind(matrix(0,nrow=nrow(beta_pena),ncol=ncol(alpha_pena)),beta_pena))


Beta_est=GetFaFQLEst(respon_Y,fact_numopt,cov_Z,pena_mat,tau,h)




# Air Pollution Analysis

respon_Y=aqi_mon
cov_W=hum_mon
cov_X=tem_dly
cov_Z=data.frame(cov_W,GetIntMat(cov_X,beta.basis))

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

fact_BIC=NULL
for(r in 1:fact_numax){
fat_num=r 
fact_res=GetFaFQLEst(respon_Y,fact_num,cov_Z,pena_mat,tau,h)
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

para_thi=para_SSE=NULL

for(j in 1:length(val_grid)){
para_thi[j,]=rep(val_grid[j],length.out=dim(X)[2])
basis_paramat=matrix(diag(para_thi[j,]),ncol=length(dim(X)[2]))
basis_penamat=matrix(basis_pena,ncol=length(basis_pena))
beta_pena=krnocker(basis_paramat,basis_penamat)
alpha_pena=matrix(0,ncol=dim(cov_W)[1],nrow=dim(cov_W)[1])
pena_mat=rbind(cbind(alpha_pena,matrix(0,nrow=nrow(alpha_pena),ncol=ncol(beta_pena))),
cbind(matrix(0,nrow=nrow(beta_pena),ncol=ncol(alpha_pena)),beta_pena))

para_res=GetFaFQLEst(respon_Y,fact_num,cov_Z,pena_mat)
para_theta=para_res$iter_theta
para_F=para_res$iter_F
para_lambda=para_res$iter_lambda

para_resYthe=para_ss=NULL
for (i in 1:N){
para_resYthe[,i]=respon_Y[,i]-cov_Z[,,i]%*%para_theta
para_ss=para_ss+sum(para_resYthe[,i]^2)
}

para_SSE[j]=para_ss
para_Fpro=martix(diag(rep(1,T)))-para_F%*%t(para_F)/T
para_Z=NULL
for(i in 1:N){
para_Z=rbind(para_Z,cov_Z[,,i])
}
para_Zpro=para_Fpro%*%para_Z
para_sthi=para_Z%*%ginv(t(para_Zpro)%*%para_Zpro+pena_mat)%*%t(para_Zpro)
para_s=martix(diag(rep(1,N*T)))-para_sthi
para_gcv[j]=para_SSE[j]/sum(diag(para_s))
}

basis_para=val_grid[which.min(para_gcv)]
basis_paramat=matrix(diag(basis_para),ncol=length(basis_para))
basis_penamat=matrix(basis_pena,ncol=length(basis_pena))
beta_pena=krnocker(basis_paramat,basis_penamat)
alpha_pena=matrix(0,ncol=dim(cov_W)[1],nrow=dim(cov_W)[1])
pena_mat=rbind(cbind(alpha_pena,matrix(0,nrow=nrow(alpha_pena),ncol=ncol(beta_pena))),
cbind(matrix(0,nrow=nrow(beta_pena),ncol=ncol(alpha_pena)),beta_pena))

Beta_est=GetFaFQLEst(respon_Y,fact_numopt,cov_Z,pena_mat,tau,h)

