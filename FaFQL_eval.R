# FaFQL model and FQLM model estimation evaluation

fafql_res=GetFaFQLEst(respon_Y,fact_num,cov_Z,pena_mat,tau,h)
fafql_theta=fafql_res$iter_theta
fafql_F=fafql_res$iter_F
fafql_lambda=fafql_res$iter_lambda

fafql_alpha=fafql_theta[1:dim(X)[2]]
fafql_gamma=fafql_theta[-(1:dim(X)[2]))]

fafqlbeta_gamma=matrix(fafql_gamma,nrow=dim(X)[2],byrow=T)
fafql_beta=fafqlbeta_gamma%*%beta_estbasis

fqlm_res=GetFQLMEst(respon_Y,cov_W,cov_B,pena_mat,tau,h)
fqlm_gamma=fqlm_res$fqlm_gamma

fqlmbeta_gamma=matrix(fqlm_gamma,nrow=dim(X)[2],byrow=T)
fqlm_beta=fqlmbeta_gamma%*%beta_estbasis

fafql_beta=fqlm_beta=fafql_srm=fqlm_srm=fafql_ssd=fqlm_ssd=NULL
fafql_Fsrm=funl_lambdasrm=NULL
for(sn in 1:sim_num){
fafql_res=GetFaFQLEst(respon_Y,fact_num,cov_Z,pena_mat,tau,h)
fafql_theta=fafql_res$iter_theta
fafql_gamma=fafql_theta[-(1:dim(X)[2]))]
fafqlbeta_gamma=matrix(fafql_gamma,nrow=dim(X)[2],byrow=T)
fafql_beta[sn,]=fafqlbeta_gamma%*%beta_estbasis

fafql_Fest=fafql_res$iter_F
fafql_lambdaest=fafql_res$iter_lambda

fafql_Fsrm[i]=apply((fafql_Fest-fafql_F)^2,2,mean)
fafql_lambdasrm[i]=apply((fafql_lambdaest-fafql_lambda)^2),2,mean)

fqlm_gamma=GetFQLMEst(respon_Y,cov_W,cov_B,pena_mat,tau,h)$fqlm_gamma
fqlmbeta_gamma=matrix(fqlm_gamma,nrow=dim(X)[2],byrow=T)
fqlm_beta[sn,]=fqlmbeta_gamma%*%beta_estbasis

fafql_srm[i]=sum((fafql_beta[i,]-beta)^2)/length(s_grid)
fqlm_srm[i]=sum((fqlm_beta[i,]-beta)^2)/length(s_grid)
}
fafql_mest=apply(fafql_beta,2,mean)
fqlm_mest=apply(fqlm_beta,2,mean)

for(s in 1:sim_num){
fafql_ssd[i]=sum((fafql_beta[i,]-fafql_mest)^2)/length(s_grid)
fqlm_ssd[i]=sum((fqlm_beta[i,]-fqlm_mest)^2)/length(s_grid)
}

fafql_rm=sqrt(mean(fafql_srm))
fqlm_rm=sqrt(mean(fqlm_srm))
fafql_sd=sqrt(mean(fafql_ssd))
fqlm_sd=sqrt(mean(fqlm_ssd))
FaFQL_Frm=sqrt(fafql_Fsrm)
FaFQL_lambdarm=sqrt(fafql_lambdasrm)




