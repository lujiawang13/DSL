Theta1procedure_f<-function(p,th){
  # Theta1
  Theta_raw_1 <- matrix(runif(p*p,min=-1, max=1), ncol=p,nrow=p) 
  Theta_raw_1<-(Theta_raw_1+t(Theta_raw_1))/2
  diag(Theta_raw_1) <- 1
  Theta_raw_1[which(Theta_raw_1>-th & Theta_raw_1<th, arr.ind = T)]<-0
  Theta1=postive_made_f(Theta_raw_1)
  return(list(Theta1=Theta1,Theta_raw_1=Theta_raw_1))
}
Theta2procedure_f<-function(Theta1procedure_result,p,sub_p,s,th){
  Theta_raw_1=Theta1procedure_result$Theta_raw_1
  indexPos<-which(Theta_raw_1>0,arr.ind = T)
  Theta_raw_1[indexPos]<-tao_PC_procedure_f(nrow(indexPos),th)
  indexNeg<-which(Theta_raw_1<0,arr.ind = T)
  Theta_raw_1[indexNeg]<-(tao_PC_procedure_f(nrow(indexNeg),th))*(-1)

  diag(Theta_raw_1) <- 1
  # Theta2
  Theta_raw_change_result<-Theta_raw_change_f(Theta_raw=Theta_raw_1,p,sub_p,s,th)
  Theta_raw_2=as.matrix(Theta_raw_change_result$Theta_class_i)
  sub_pos=Theta_raw_change_result$sub_pos
  Theta2=postive_made_f(Theta_raw_2)
  return(list(Theta2=Theta2,Theta_raw_2=Theta_raw_2,sub_pos=sub_pos))
}
Xprocedure_f<-function(Theta,p,ni,N,h){
  #################################################################
  mu=rep(0,p)
  I=diag(rep(0,p));dim(I)
  #########################################################
  Xp=c();
  xList<-list()
  for(i in 1:N){
    Thetaitemp=rmvnorm(h,mean = mu, sigma = Theta, method=c("eigen", "svd", "chol"), pre0.9_9994 = FALSE)
    Thetai=t(Thetaitemp)%*%Thetaitemp;

    Sigmai=solve(Thetai)
    
    Xi=rmvnorm(ni, mean = mu, sigma = Sigmai, method=c("eigen", "svd", "chol"), pre0.9_9994 = FALSE)
    Xp=rbind(Xp,Xi)
    xList[[i]]<-Xi
  }

  return(list(Xp=Xp,xList=xList))
}



