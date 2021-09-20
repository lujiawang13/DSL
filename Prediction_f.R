Prediction_f<-function(pos,TestSi_List,Testnivector,AllTheta1,AllTheta2,p,h1,h2){
  if(identical(as.character(pos),character(0))==FALSE){
    ##############sub  Sigma test
    p_sub=length(pos)
    sub_TestSigma_List<-list()
    for(ind in 1:length(Testnivector)){
      Si=TestSi_List[[ind]]
      sub_Si<-Si[pos,pos]
      sub_TestSigma_List[[ind]]<-sub_Si
    }
    ###################Shared Theta of class 1 and class 2
    sub_AllTheta1<-AllTheta1[pos,pos]
    sub_AllTheta2<-AllTheta2[pos,pos] 
    
    #################Predition for test
    num=length(Testnivector);nnewvector=Testnivector
    Prediction_output<-Prediction_AUC(sub_TestSigma_List,nnewvector,sub_AllTheta1, sub_AllTheta2,p_sub,h1,h2)
    predprob<-Prediction_output$prob_P_AUC_chain
    return(predprob)
  }else{
    return("Error in pos")
  }
}


################################Prediction 
#num: number of dataSi
Prediction_AUC<-function(dataSi_List,nnewvector,AllTheta1, AllTheta2,p,h1,h2){
  # dataSigma=sub_Sigma1;p=p_sub;AllTheta1=sub_AllTheta1;AllTheta2= sub_AllTheta2;h1=sub_h_Left;h2=sub_h2
  # dataSigma=sub_Sigma_test;p=p_sub;AllTheta1=sub_AllTheta1;AllTheta2= sub_AllTheta2;h1=sub_h_left;h2=sub_h2
  num=length(nnewvector)
  
  # class 1
  AllTheta1=as.matrix(AllTheta1); 
  detbP1=(-h1/2)*log(det(AllTheta1))
  # if(is.singular.matrix(AllTheta1, tol = 1e-08)){AllTheta1=diag(1,p)}
  # class 2
  AllTheta2=as.matrix(AllTheta2);
  detbNC1=(-h2/2)*log(det(AllTheta2))
  # if(is.singular.matrix(AllTheta2, tol = 1e-08)){AllTheta2=diag(1,p)}
  
  determineP_log=c();
  prob_P_AUC_chain<-c()
  for(i in 1:num){        
    Si=dataSi_List[[i]];  
    Si=as.matrix(Si) 
    nnew<-nnewvector[i]  
    
    ######log
    # class 1
    detbP2=(-(nnew+h1-1)/2)*log(det(nnew*as.matrix(Si)+solve(AllTheta1)));
    # class 2
    detbNC2=(-(nnew+h2-1)/2)*log(det(nnew*as.matrix(Si)+solve(AllTheta2)));
    
    GammaP_P=0;GammaP_NC=0;
    for (k in 1:p){
      # class 1
      gammab=lgamma((nnew+h1-k)/2)-lgamma((h1-k+1)/2);
      GammaP_P=GammaP_P+gammab;
      # class 2
      gammaNCb=lgamma((nnew+h2-k)/2)-lgamma((h2-k+1)/2);
      GammaP_NC=GammaP_NC+gammaNCb;
    }	
    
    log_prob_P_AUC<-detbP1+detbP2+GammaP_P
    log_prob_NC_AUC<-detbNC1+detbNC2+GammaP_NC
    #prob_P_AUC=log_prob_NC_AUC/(log_prob_P_AUC+log_prob_NC_AUC)  
    prob_P_AUC=1/(1+exp(log_prob_NC_AUC-log_prob_P_AUC))
    prob_P_AUC_chain<-c(prob_P_AUC_chain,prob_P_AUC)
    
    Comp_log=log_prob_P_AUC-log_prob_NC_AUC
    determineP_log=c(determineP_log,Comp_log)
  }
  
  result_with_log=sum(determineP_log>0)
  
  return(list(result_with_log=result_with_log,Seqof_Prob_log=determineP_log,
              prob_P_AUC_chain=prob_P_AUC_chain))
  
  
}