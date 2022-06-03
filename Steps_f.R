###tao produce
tao_PC_procedure_f<-function(Num,th){
  tao=runif(Num, min=th, max=1)
  return(tao)
}

###tao produce
tao_procedure_f<-function(Num,th){
 tao_temp=runif(Num, min=-th, max=th)
 tao=tao_temp
 tao[which(tao_temp<0)]=tao[which(tao_temp<0)]/th*(1-th)-th
 tao[which(tao_temp>0)]=tao[which(tao_temp>0)]/th*(1-th)+th
 return(tao)
}
taoRe_procedure_f<-function(value,th){
  # delta=0.01
  delta=0.05
  if(value<0){
    tao=runif(1, min=max(-1,value-delta), max=min(-0.5,value+delta))
  }
  if(value>=0){
    tao=runif(1, min=max(0.5,value-delta), max=min(1,value+delta))
  }

  return(tao)
}

#####vectorize_upper_triangle
vectorize_upper_triangle_f<-function(Theta){
  p=ncol(Theta)
  vector_chain<-c()
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      vector_chain<-c(vector_chain,Theta[i,j])
    }
  }
  return(vector_chain)
}
######change to marix
recover_to_matrix_f<-function(original_Theta,vector_chain){
  p=ncol(original_Theta)
  recovermatrix<-matrix(NA,ncol=p,nrow=p)
  k=0
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      k=k+1
      recovermatrix[i,j]=vector_chain[k]
      recovermatrix[j,i]=vector_chain[k]
    }
  }
  for(i in 1:p){
    recovermatrix[i,i]=original_Theta[i,i]
  }
  return(recovermatrix)
}


##### Step 2.3
Theta_raw_change_f<-function(Theta_raw,p,sub_p,s,th){
  sub_pos<-sort(sample(1:p,sub_p,replace = F))
  
  sub_Theta<-Theta_raw[sub_pos,sub_pos]
  sub_Theta_vector<-vectorize_upper_triangle_f(Theta=sub_Theta)
  
  ####Step 2
  nonZero_num=sum(sub_Theta_vector!=0)
  YesZero_num=sum(sub_Theta_vector==0)
  NumModify=ceiling(s*nonZero_num)
  if(YesZero_num<NumModify){
    NumModify=YesZero_num
    print(paste("nonZero_num",nonZero_num,"YesZero_num",YesZero_num))
  }
  
  #########Non-zero to zero
  Non_Zero_pos=which(sub_Theta_vector!=0)
  if(length(Non_Zero_pos)==1){
    Non_Zero_pos_select=Non_Zero_pos
  }else{
    Non_Zero_pos_select=sample(Non_Zero_pos,NumModify,replace=F)
  }
  
  ### zero to non-zero
  Yes_Zero_pos=which(sub_Theta_vector==0)
  if(length(Yes_Zero_pos)==1){
    Yes_Zero_pos_select=Yes_Zero_pos
  }else{
    Yes_Zero_pos_select=sample(Yes_Zero_pos,NumModify,replace=F)
  }
  Newvalue=tao_procedure_f(Num=NumModify,th)
  
  ###Step 2.4 for the remaining
  Theta_vector<-vectorize_upper_triangle_f(Theta=Theta_raw)
  Non_Zero_posAll=which(Theta_vector!=0)
  # Newvalue_all=tao_procedure_f(Num=length(Non_Zero_posAll),th)
  Newvalue_all=sapply(Theta_vector[Non_Zero_posAll],taoRe_procedure_f,th)
  
  sub_Theta_vector[Non_Zero_posAll]=Newvalue_all;
  sub_Theta_vector[Non_Zero_pos_select]=0;
  sub_Theta_vector[Yes_Zero_pos_select]=Newvalue;
  # 
  
  ######change to marix
  sub_Theta_class_i<-recover_to_matrix_f(original_Theta=sub_Theta,vector_chain=sub_Theta_vector)
  Theta_raw[sub_pos,sub_pos]=sub_Theta_class_i
  Theta_class_i=Theta_raw
  return(list(Theta_class_i=Theta_class_i,sub_pos=sub_pos))
}


# Theta_raw_change_f<-function(Theta_raw,p,sub_p,s,th){
#   sub_pos<-sort(sample(1:p,sub_p,replace = F))
#   sub_Theta<-Theta_raw[sub_pos,sub_pos]
#   
#   ####Step 2
#   nonZero_num=sum(sub_Theta!=0)
#   YesZero_num=sum(sub_Theta==0)
#   NumModify=ceiling(s*nonZero_num)
#   if(YesZero_num<NumModify){
#     NumModify=YesZero_num
#     print(paste("nonZero_num",nonZero_num,"YesZero_num",YesZero_num))
#   }
#   
#   #########Non-zero to zero
#   Non_Zero_pos=which(sub_Theta!=0)
#   if(length(Non_Zero_pos)==1){
#     Non_Zero_pos_select=Non_Zero_pos
#   }else{
#     Non_Zero_pos_select=sample(Non_Zero_pos,NumModify,replace=F)
#   }
#   
#   ### zero to non-zero
#   Yes_Zero_pos=which(sub_Theta==0)
#   if(length(Yes_Zero_pos)==1){
#     Yes_Zero_pos_select=Yes_Zero_pos
#   }else{
#     Yes_Zero_pos_select=sample(Yes_Zero_pos,NumModify,replace=F)
#   }
#   Newvalue=tao_procedure_f(Num=NumModify,th)
#   
#   sub_Theta[Non_Zero_pos_select]=0;
#   sub_Theta[Yes_Zero_pos_select]=Newvalue;
#   
#   
#   ###Step 2.4 for the remaining
#   Non_Zero_posAll=which(Theta_raw!=0)
#   Newvalue_All=tao_procedure_f(Num=length(Non_Zero_posAll),th)
#   
#   Theta_raw[Non_Zero_posAll]=Newvalue_All;
#   Theta_raw[sub_pos,sub_pos]=sub_Theta
#   Theta_class_i=Theta_raw
#   return(list(Theta_class_i=Theta_class_i,sub_pos=sub_pos))
# }

### Step 3
postive_made_f<-function(Theta){
  p=ncol(Theta)
  #############ensure positive
  Thetaabsolute<-abs(Theta);Thetaabsolute
  ThetaabsoluteRow=rowSums(Thetaabsolute)-1
  Div=1.5*ThetaabsoluteRow
  for(i in 1:p){
    for(j in 1:p){
      if(i!=j){
        Theta[i,j]=Theta[i,j]/Div[i]}
    }
  }
  Theta_positive=(t(Theta)+Theta)/2
  return(Theta_positive)
}


