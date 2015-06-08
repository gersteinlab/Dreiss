options(stringsAsFactors=FALSE)
eigendecomp=function(data,r=5,flag_msg=T){
  svd1=svd(data);U1=svd1$u;S1=svd1$d;V1=svd1$v
  if(flag_msg){svpecent=(S1^2/sum(S1^2))[1:r];print(svpecent);print(sum(svpecent))}
  return(list(W=U1[,1:r]%*%diag(S1[1:r]),E=t(V1[,1:r])))
}
eigensysmat=function(X,U,r=5,r2=r,flag_msg=T){
  Xeigen=eigendecomp(X,r,flag_msg);Wx=Xeigen$W;Xr=Xeigen$E
  Ueigen=eigendecomp(U,r2,flag_msg);Wu=Ueigen$W;Ur=Ueigen$E
  rownames(Wx)=rownames(X);rownames(Wu)=rownames(U)
  X2=Xr[,2:ncol(X)];Y2=rbind(Xr[,1:(ncol(Xr)-1)],Ur[,1:(ncol(Ur)-1)])
  Z=X2%*%pseudoinverse(Y2)
  A=Z[,1:r];B=Z[,(r+1):ncol(Z)]   # A and B of eigensystem
  Aeigen=eigen(A);Beigen=eigen(B) #eigenvalues/vectors of A and B
  Aev=Aeigen$values;Bev=Beigen$values
  Aevec=Aeigen$vectors;Bevec=Beigen$vectors
  Acoeff=Wx%*%Aevec;rownames(Acoeff)=rownames(X) #gene coeffs over principal dynamic patterns on orthologous space
  Bcoeff=Wx%*%Bevec;rownames(Bcoeff)=rownames(X) #gene coeffs over principal dynamic patterns on regulation space  
  return(list(A=A,B=B,Aev=Aev,Bev=Bev,Aevec=Aevec,Bevec=Bevec,
              Acoeff=Acoeff,Bcoeff=Bcoeff,Xr=Xr,Ur=Ur,Wx=Wx,Wu=Wu))
}