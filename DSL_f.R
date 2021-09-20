####Si
dataSi_f<-function(X_List,nivector,p){
  n=length(X_List)
  dataSi_List<-list();
  for(i in 1:n){
    ni=nivector[i]
    datai=X_List[[i]];      
    bar_datai=array(colMeans(datai),c(p,1));
    Ai=t(datai)%*%as.matrix(datai)-ni*bar_datai%*%t(bar_datai);
    Si=Ai/ni;
    dataSi_List[[i]]<-Si
  }
  bar=bar_datai%*%t(bar_datai)
  return(dataSi_List);
}

# Train
Theta_f<-function(TrainSi_List,Trainnivector,TrainY,p,h1,h2,lamda1,lamda2,mu1,mu2,area_k,
                       Opt_method,structure,
                  betamethod){
  # initial
  ntimes=1
  AOiterUsed=20
  N1=sum(TrainY==1)
  N2=sum(TrainY==0)
  lamda1=lamda1/N1
  lamda2=lamda2/N2
  mu1=mu1/N1
  mu2=mu2/N2
  # AO
  Index_Left<-which(TrainY==1)
  Index_Right<-which(TrainY==0)
  AllTheta1<-InitialTheta_f(dataSi_List=TrainSi_List[Index_Left],nivector=Trainnivector[Index_Left],p,h1)
  AllTheta2<-InitialTheta_f(dataSi_List=TrainSi_List[Index_Right],nivector=Trainnivector[Index_Right],p,h2)
  AllTheta10=diag(p);AllTheta20=diag(p)
  pos_c<-rep(1,p)
  
  #######################
  best_pos=1:p

  
  ###############################
  AOiter=0
  while(sum(abs(AllTheta10-AllTheta1)+abs(AllTheta20-AllTheta2))>0.001){
    AOiter=AOiter+1
    print(paste("AOiter",AOiter))
    AllTheta10=AllTheta1;AllTheta20=AllTheta2
    TrainResult<-OptAO_f(TrainSi_List,Trainnivector,TrainY,AllTheta1,AllTheta2,pos_c,p,h1,h2,
                         lamda1,lamda2,mu1,mu2,area_k,
                         Opt_method,ntimes,structure,
                         betamethod)
    AllTheta1=TrainResult$AllTheta1
    AllTheta2=TrainResult$AllTheta2
    best_pos=TrainResult$best_pos
    pos_c=rep(0,p)
    pos_c[best_pos]=1
    
    if(AOiter==AOiterUsed) break
    
   
  }  
  return(list(AllTheta1=AllTheta1,AllTheta2=AllTheta2,
              best_pos=best_pos,pos_c=pos_c))
}
loggammap_f<-function(n,p){
  A=n/2-((1:p)-1)/2
  loggammap=p*(p-1)/4*log(pi)+sum(log(gamma(A)))
  # gammap=exp(loggammap)
  return(loggammap)
}

OptAO_f<-function(TrainSi_List,Trainnivector,TrainY,AllTheta1,AllTheta2,pos_c,p,h1,h2,
                  lamda1,lamda2,mu1,mu2,area_k,
                  Opt_method,ntimes,structure,
                  betamethod){
  Index_Left<-which(TrainY==1)
  Index_Right<-which(TrainY==0)
  
  AllTheta1<-ThetaEM_f(dataSi_List=TrainSi_List[Index_Left],nivector=Trainnivector[Index_Left],Theta0=AllTheta1,Z=AllTheta2,pos_c,
                           p,h1,lamda1,mu1,
                           betamethod)
  
  AllTheta2<-ThetaEM_f(dataSi_List=TrainSi_List[Index_Right],nivector=Trainnivector[Index_Right],Theta0=AllTheta2,Z=AllTheta1,pos_c,
                            p,h2,lamda2,mu2,
                            betamethod)
  
  pos_result<-SubGraph_f(AllTheta1,AllTheta2,p,area_k,Opt_method,ntimes,structure)
  pos_matrix=pos_result$pos_matrix
  best_pos=pos_result$best_pos
  if(is.null(dim(pos_matrix))){pos_matrix<-matrix(pos_matrix,nrow=1)}

  
  return(list(AllTheta1=AllTheta1,AllTheta2=AllTheta2,
              pos_matrix=pos_matrix,best_pos=best_pos))
}

