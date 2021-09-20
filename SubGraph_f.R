SubGraph_f<-function(AllTheta1,AllTheta2,p,area_k,Opt_method,ntimes=NA,structure){
  ##### binray structure difference
  if(structure=='binarystructure'){
    Theta_Stucture_left<-matrix(as.numeric(AllTheta1!=0),c(p,p))
    Theta_Stucture_right<-matrix(as.numeric(AllTheta2!=0),c(p,p))
    matrix_diff<-abs(Theta_Stucture_left-Theta_Stucture_right)
  }
  # weighted structure difference
  if(structure=='weightedstructure'){
    matrix_diff<-abs(AllTheta1-AllTheta2)
  }
  ##########sub pos
  pos_result<-Opt_subPos_f(matrix_diff,area_k,p,ntimes,Opt_method)
  pos_matrix=pos_result$pos_matrix 
  best_pos=t(as.matrix(pos_result$best_pos))
  
  return(list(pos_matrix=pos_matrix,best_pos=best_pos))
}

Opt_subPos_f<-function(matrix_diff,area_k,p,ntimes,Opt_method){
  if(Opt_method=="Rsymphony"){
    if(is.na(ntimes)){ntimes=50}
    pos_matrix<-c()
    sol_vector_sum=0
    for(i in 1:ntimes){
      opt<-sub_Extract_Rsymphony_pos(matrix_diff,area_k,p)
      sol_vector<-opt$sol_vector
      sol_vector_sum=sol_vector_sum+sol_vector
      pos_temp=which(sol_vector==1); 
      Check_Repeat<-Check_Repeat_f(pos_matrix,pos_temp)
      if(Check_Repeat==TRUE){
        next
      }
      pos_matrix<-rbind(pos_matrix,pos_temp)
    }
  }else if(Opt_method=="RCplex"){
    if(Constraints=="Constraints"){opt_result<-RCplexConstraints_f(matrix_diff,area_k,p,ntimes)}else{
      opt_result<-RCplex_f(matrix_diff,area_k,p,ntimes)}
    pos_matrix<-c()
    sol_vector_sum=0
    for(i in 1:length(opt_result)){
      sol<-opt_result[[i]]
      sol_vector<-sol$xopt #obj<-sol$obj
      sol_vector<-sol_vector[1:p]
      sol_vector_sum=sol_vector_sum+sol_vector
      pos_temp=which(sol_vector>0.5); 
      if(length(opt_result)>1 & length(pos_temp)!=area_k){next}
      pos_matrix<-rbind(pos_matrix,pos_temp)
    }
  }
  best_pos<-which(sol_vector_sum >= sort(sol_vector_sum, decreasing=T)[area_k])
  if(length(best_pos)>area_k){best_pos=best_pos[1:area_k]}
  best_pos=pos_matrix[sample(1:nrow(pos_matrix),1),]
  return(list(pos_matrix=pos_matrix,best_pos=best_pos))
}




Check_Repeat_f<-function(pos_matrix,pos_temp){
  Check_Repeat=FALSE
  if(is.null(pos_matrix)==FALSE){
    for(rr in 1:nrow(pos_matrix)){
      if(all(pos_matrix[rr,]==pos_temp)){
        Check_Repeat=TRUE
        break
      }
    }
  }
  return(Check_Repeat)
  
}

sub_Extract_Rsymphony_pos<-function(matrix_diff,area_k,p){
  
  ###Constraint (p+p*p)*(sum(matrix_diff)/2+1)
  # num_1=sum(matrix_diff)/2
  ep=1e-6
  num_1=0
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      if(abs(matrix_diff[i,j])>ep){
        num_1=num_1+1
      }
    }
  }
  if(floor(num_1)!=num_1){
    print("Error not integer in structures!")
    # break;
  }
  A.con_1<-matrix(0,ncol=p+num_1,nrow=num_1)
  obj=c()
  pos=0
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      if(abs(matrix_diff[i,j])>ep){
        pos=pos+1
        A.con_1[pos,i]=1
        A.con_1[pos,j]=1
        A.con_1[pos,p+pos]=-2
        obj=c(obj,matrix_diff[i,j])
      }
    }
  }
  
  direction ="max"
  objective.in=c(diag(matrix_diff),2*obj)
  A.con_0=-A.con_1
  const.mat_temp=rbind(A.con_1,A.con_0)
  const.mat=rbind(const.mat_temp,c(rep(1,p),rep(0,num_1)))
  const.rhs<-c(rep(1,num_1),rep(0,num_1),area_k)
  const.dir<-rep("<=",length(const.rhs))
  varible.type<-rep("B",p+num_1)
  sol=Rsymphony_solve_LP(obj=objective.in, mat=const.mat, dir=const.dir, time_limit =3*60,node_limit = 300,
                         rhs=const.rhs,types =varible.type, max = TRUE)
  
  
  status=sol$status
  sol_vector<-sol$solution
  ##############
  #pos_temp=which(sol_vector==1)[1:18]; pos_temp
  
  # install.packages("lpSolve")
  # library(lpSolve)
  # sol=lp (direction = "max", objective.in, const.mat, const.dir, const.rhs,
  #         transpose.constraints = TRUE, presolve=0, compute.sens=0,
  #         binary.vec=1:(p+num_1), all.int=FALSE, all.bin=FALSE, scale = 196,
  #         num.bin.solns=1, use.rw=FALSE)
  
  
  # if(status>0){
  #   print("No feasible solution")
  #   break;
  # }
  obj<-sol$objval
  return(list(obj=obj,sol_vector=sol_vector[1:p],status=status))
  
}

#########################################################
RCplex_f<-function(matrix_diff,area_k,p,ntimes){
  Qmat=matrix_diff
  cvec=matrix(rep(0,p),nrow=p,ncol=1)
  Amat <- matrix(rep(1,p),nrow=1,ncol=p)
  bvec <- area_k+0.0
  vtype <- rep("B",p)
  sol <- Rcplex(cvec,Amat,bvec,Qmat,objsense="max",sense=c('L'),vtype=vtype,n=ntimes);
  # status=sol$status
  # sol_vector<-sol$xopt
  # obj<-sol$obj
  return(sol)
  
}


############################
# -x1-x2<=-1.5+M*y
# x1+x2<=0.5+M*(1-y)

###
# -x1-x2-M*y<=-1.5
# x1+x2+M*y<=0.5+M
RCplexConstraints_f<-function(matrix_diff,area_k,p,ntimes){
  Equ_matrix<-rbind(t(matrix(1:30,nrow=2,ncol=15)),c(32,33))
  nforEqu=nrow(Equ_matrix)
  
  Qmat=matrix_diff
  cvec=matrix(rep(0,p+nforEqu),nrow=(p+nforEqu),ncol=1)
  Amat <- matrix(c(rep(1,p),rep(0,nforEqu)),nrow=1,ncol=(p+nforEqu))
  bvec <- area_k+0.0
  vtype <- rep("B",(p+nforEqu))
  ###########################################
  M=10
  bvecforEqu<-c(rep(-1.5,nforEqu),rep(0.5+M,nforEqu))
  AmatforEqu1<-matrix(0,nrow=nforEqu,ncol=(p+nforEqu))
  AmatforEqu2<-matrix(0,nrow=nforEqu,ncol=(p+nforEqu))
  for(i in 1:nforEqu){
    pos=Equ_matrix[i,]
    AmatforEqu1[i,pos]=-1
    AmatforEqu1[i,(p+i)]=-M
    AmatforEqu2[i,pos]=1
    AmatforEqu2[i,(p+i)]=M
  }
  Amat<-rbind(Amat,AmatforEqu1,AmatforEqu2)
  bvec<-c(bvec,bvecforEqu)
  Qmat<-rbind(Qmat,matrix(0,nrow=nforEqu,ncol=p))
  Qmat<-cbind(Qmat,matrix(0,nrow=(p+nforEqu),ncol=nforEqu))
  sol <- Rcplex(cvec,Amat,bvec,Qmat,objsense="max",sense=c('L'),vtype=vtype,n=ntimes);
  return(sol)
  
}









