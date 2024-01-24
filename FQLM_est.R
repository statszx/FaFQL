# calculate the FQLM estimates

GetFQLMEst=function(respon_Y,cov_W,cov_B,pena_mat,tau,h){
cov_Z=covZ_mat=covZY_mat=NULL
for (i in 1:N){
cov_Z[,,i]=cbind(cov_W[,,i],cov_B[,,i])
covZ_mat=covZ_mat+t(cov_Z[,,i])%*%cov_Z[,,i]
covZY_mat=covZY_mat+t(cov_Z[,,i])%*%respon_Y[,i]
}

resi_Y=NULL
smorho=NULL
covZ_mat=covZY_mat=NULL
for (i in 1:N){
resi_Y[,i]=respon_Y[,i]
covZ_mat=covZ_mat+t(cov_Z[,,i])%*%cov_Z[,,i]
covZY_mat=covZY_mat+t(cov_Z[,,i])%*%resi_Y[,i]
smorho=smorho+GetSmoRho(resi_Y[,i],tau,h)
}

fqlm_res=smorho+t(iter_theta)%*%pena_mat%*%iter_theta*pena_para)

fqlm_ge=fqlm_e=NULL
for (i in 1:n){
fqlm_ge=fqlm_ge+GetGauKer(respon_Y[i]-sum(cov_Z[i,]*iter_theta),h)*
as.matrix(c(1,cov_Z[i,]))%*%c(1,cov_Z[i,])
fqlm_e=fqlm_e-(GauKer(respon_Y[i]-sum(cov_Z[i,]*iter_theta),h)+tau-1)*
c(1,cov_Z[i,])}
fqlm_e=as.matrix(fqlm_e)
fqlm_matrix=bdiag(matrix(0,nrow=nrow(cov_Z),ncol=ncol(cov_Z)+1),
pena_para*pena_mat)+fqlm_ge
fqlm_mat=svd(fqlm_matrix)
fqlm_matv=(fqlm_mat$u)%*%diag((fqlm_mat$d)^-1)%*%t(fqlm_mat$v)
fqlm_matve=fqlm_matv%*%fqlm_e
fqlm_theta=fqlm_matve

fqlm_alpha=fqlm_theta[1:dim(X)[2]]
fqlm_gamma=fqlm_theta[-(1:dim(X)[2]))]

beta_gamma=matrix(fqlm_gamma,nrow=dim(X)[2],byrow=T)
fqlm_beta=beta_gamma%*%beta_estbasis

out=list(fqlm_beta,fqlm_gamma,fqlm_theta)
return(out)
}


