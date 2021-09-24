ThetaEM_f<-function(dataSi_List,nivector,Theta0,Z,pos_c,p,h,lamda,mu,
                    betamethod){
  ######share Theta 
  ThetaEM=EM_f(dataSi_List,nivector,Theta0,Z,pos_c,p,h,lamda,mu,
                        betamethod)
  AllTheta=ThetaEM$Thetat

  return(AllTheta)
}


###Initialization value
InitialTheta_f<-function(dataSi_List,nivector,p,h){
  n=length(nivector)
  Ai_sum<-matrix(0,nrow=p,ncol=p)
  for(mm in 1:length(nivector)){ 
    Si=dataSi_List[[mm]];
    Ai=solve(Si)
    Ai_sum=Ai_sum+Ai*(1-1/nivector[mm])
  }
  Theta0=Ai_sum/n/h
  return(Theta0)
}
  

EM_f<-function(dataSi_List,nivector,Theta0,Z,pos_c,p,h,lamda,mu,
                        betamethod){
  n=length(nivector)
  Theta2t=Theta0
  #######
  # optEM_chain<-c()
  j=0;
  repeat{
    j=j+1;
    Spie=0;
    for(i in 1:n){
      ni<-nivector[i];
      Si=dataSi_List[[i]];
      Di=(ni+h-1)*solve((ni*Si+solve(Theta2t)));
      Spie=Spie+Di/n/h;
    }
    
  
    ###desent
    BCD_result<-BCD_f(H=Spie,Theta0=Theta2t,Z,pos_c,lamda,mu,epsilon=1e-6,
                      betamethod)
    Theta2t0 <- BCD_result$Theta;
    
    if (abs(norm_f(Theta2t0-Theta2t,'2'))<1e-8 | j>50){break;}
    Theta2t=Theta2t0; 
  }
  
  return(list(Thetat = Theta2t0, Spie = Spie, iter = j))
}

norm_f<-function(Theta,type){
  if(type=='1'){
    v=sum(abs(Theta))
  }
  if(type=='2'){
    v=sum(Theta^2)
  }
  return(v)
}
