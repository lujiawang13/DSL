######   Solve Discriminant Subgraph Learning (DSL) (Wang, et al., 2021)
#   Training Input:
#   TrainX_List: Training data X
#   TrainY: binary Y
#   h1/h2: parameters of the Wishart distribution in Class 1/class 2
#   lamda1/lamda2: tuning paramters
#   mu1/mu2: tuning paramters
#   area_k: estimated column/row number of the discriminant subgraph
#   Opt_method: method of solving quadratic programming. Opt_method="RCplex" or "Rsymphony" 
#   structure: weighted IC difference or structure difference; structure='binarystructure' or 'weightedstructure'
#   betamethod: method to estimate one column/row in Theta;betamethod='majorizationminimization';betamethod='cd';betamethod='others'

#   Training Output:
#   EstimatedTheta1: estimated graph in class 1
#   EstimatedTheta2: estimated graph in class 2
#   best_pos: estimated subgraph position

#   Testing Input:
#   TestX_List: Testing data X
#   Output:
#   predprob: predictive probability to be class 1

# Written by Lujia Wang @ Georgia Tech


# library(kernlab)
# library(Rsymphony)
# library(Rcplex)
# library(AUC)
#############################################
rm(list=ls(all=TRUE));
root="C:\\" # path
source(paste(root,"BCD_f.R",sep=""))
source(paste(root,"BCD_para_f.R",sep=""))
source(paste(root,"EM_f.R",sep=""))
source(paste(root,"DSL_f.R",sep=""))
source(paste(root,"SubGraph_f.R",sep=""))
source(paste(root,"Prediction_f.R",sep=""))
############################### data reading

# Training Input
Datapath0<-paste(root,"Simulation data\\",sep='')
Datapath<-paste(Datapath0,"Train data\\",sep='')

datafile<-paste(Datapath,'Ybinary.csv',sep='')
TrainY<-read.csv(datafile)[,2]

TrainX_List<-list()
Trainnivector<-c()
for(i in 1:length(list.files(paste(Datapath,'Xdata\\',sep='')))){
  datafile<-paste(Datapath,'Xdata\\Sample_',i,'.csv',sep='')
  data<-read.csv(datafile)
  TrainX_List[[i]]=data
  Trainnivector<-c(Trainnivector,nrow(data))
}
p=ncol(data);p


#####################################################
# data preproses
TrainSi_List<-dataSi_f(X_List=TrainX_List,nivector=Trainnivector,p)


# DSL training
Opt_method="RCplex"
# Opt_method="Rsymphony" 
Constraints="NoConstraints"
# structure='binarystructure'
structure='weightedstructure'
# betamethod='majorizationminimization'
betamethod='cd'

h=100;
lamda=2.5;
mu=0.01
h1=h;h2=h
lamda1=lamda;lamda2=lamda
mu1=mu;mu2=mu
area_k=40

###############  
# Train
Theta_Info<-Theta_f(TrainSi_List=TrainSi_List,Trainnivector=Trainnivector,TrainY=TrainY,p,h1,h2,lamda1,lamda2,mu1,mu2,area_k,
                    Opt_method=Opt_method,structure=structure,
                    betamethod=betamethod)
EstimatedTheta1=Theta_Info$AllTheta1
EstimatedTheta2=Theta_Info$AllTheta2
best_pos=Theta_Info$best_pos  

# structural comparison
# True sub-graph positions and True Theta1 & True Theta2
True_sub_pos<-read.csv(paste(Datapath,'sub_pos.csv',sep=''))[,1]
True_Theta1<-read.csv(paste(Datapath,'Theta1.csv',sep=''))
True_Theta2<-read.csv(paste(Datapath,'Theta2.csv',sep=''))

# comparison
estimate_pos=best_pos
intersect_pos<-intersect(estimate_pos,True_sub_pos)
length(intersect_pos) # accurately identified subgraph 


# Structrual accuracy for True Theta1 vs Estimated Theta1
EstimatedTheta1_Structure<-matrix(as.numeric(EstimatedTheta1!=0),c(p,p))
TrueTheta1_Structure=matrix(as.numeric(True_Theta1!=0),c(p,p))
Truevector1<-TrueTheta1_Structure[lower.tri(TrueTheta1_Structure,diag = FALSE)]
Estimatedvector1<-EstimatedTheta1_Structure[lower.tri(EstimatedTheta1_Structure,diag = FALSE)]
StructuralAcc1<-sum(Truevector1==Estimatedvector1)/(p*(p-1)/2)

# Structrual accuracy for True Theta2 vs Estimated Theta2
EstimatedTheta2_Structure<-matrix(as.numeric(EstimatedTheta2!=0),c(p,p))
TrueTheta2_Structure=matrix(as.numeric(True_Theta2!=0),c(p,p))
Truevector2<-TrueTheta2_Structure[lower.tri(TrueTheta2_Structure,diag = FALSE)]
Estimatedvector2<-EstimatedTheta2_Structure[lower.tri(EstimatedTheta2_Structure,diag = FALSE)]
StructuralAcc2<-sum(Truevector2==Estimatedvector2)/(p*(p-1)/2)

StructuralAccAve=(StructuralAcc1+StructuralAcc2)/2;StructuralAccAve


############### Test Input
# Test data
TestX_List<-list()
Testnivector<-c()
Datapath<-paste(Datapath0,"Test data\\",sep='')
for(i in 1:length(list.files(paste(Datapath,'Xdata\\',sep='')))){
  datafile<-paste(Datapath,'Xdata\\TestSample_',i,'.csv',sep='')
  data<-read.csv(datafile)
  TestX_List[[i]]=data
  Testnivector<-c(Testnivector,nrow(data))
}
# data preproses
TestSi_List<-dataSi_f(X_List=TestX_List,nivector=Testnivector,p)

predprob<-Prediction_f(best_pos,TestSi_List=TestSi_List,Testnivector=Testnivector,
                            AllTheta1=EstimatedTheta1,AllTheta2=EstimatedTheta2,p,h1,h2)

# Test Accuracy  
TestY=read.csv(paste(Datapath,'TestYbinary.csv',sep=''))[,'Y']

overall_predict_AUC=auc(roc(predprob,factor(TestY)));overall_predict_AUC

