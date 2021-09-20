######   Solve Discriminant Subgraph Learning (DSL) (Wang, et al., 2021)
#   Input:
#   TrainX_List: data X
#   TrainY: binary Y
#   Trainnivector: a chain of ni to each Xi in TrainX_List
#   X1_List: data list in Class 1
#   X2_List: data list in Class 2
#   h1/h2: parameters of the Wishart distribution in Class 1/class 2
#   lamda1/lamda2: tuning paramters
#   mu1/mu2: tuning paramters
#   area_k: estimated column/row number of the discriminant subgraph
#   Opt_method: method of solving quadratic programming
#   structure: weighted graph difference or structure difference.structure='binarystructure' or 'weightedstructure'
#   betamethod: method to estimate one column/row in Theta

#   Output:
#   EstimatedTheta1: estimated graph in class 1
#   EstimatedTheta2: estimated graph in class 2
#   best_pos: estimated subgraph position
#   predprob: predictive probability to be class 1

# Written by Lujia Wang @ Georgia Tech


# library(kernlab)
# library(Rsymphony)
# library(Rcplex)
#############################################
rm(list=ls(all=TRUE));
root="C:\\Users\\lwang724\\Dropbox (GaTech)\\IC code\\Intergration\\github\\"
source(paste(root,"Steps_f.R",sep=''))
source(paste(root,"DataProcedure_f.R",sep=''))
source(paste(root,"BCD_f.R",sep=""))
source(paste(root,"BCD_para_f.R",sep=""))
source(paste(root,"EM_f.R",sep=""))
source(paste(root,"DSL_f.R",sep=""))
source(paste(root,"SubGraph_f.R",sep=""))
source(paste(root,"Prediction_f.R",sep=""))
############################### simulation data generation (See 'Data produce Run.R' in detail)
ni=200
N1=50
N1=100
N2=N1
N1test=100
N2test=100
p=50;sub_p=40
h1=100
h2=h1
s=0.5
th=0.5

Theta1procedure_result<-Theta1procedure_f(p,th)
Theta1=Theta1procedure_result$Theta1


Theta2procedure_result<-Theta2procedure_f(Theta1procedure_result,p,sub_p,s,th)
Theta2=Theta2procedure_result$Theta2

True_sub_pos=Theta2procedure_result$sub_pos

# Train and Test data
data1_result<-Xprocedure_f(Theta=Theta1,p,ni,N=N1+N1test,h=h1)
X1_List<-data1_result$xList

data2_result<-Xprocedure_f(Theta=Theta2,p,ni,N=N2+N2test,h=h2)
X2_List<-data2_result$xList

# Train data
TrainX_List<-c(X1_List[1:N1],X2_List[1:N2])
TrainY<-c(rep(1,N1),rep(0,N2))
Trainnivector<-rep(ni,(N1+N2))

# Test data
TestX_List<-c(X1_List[(N1+1):(N1+N1test)],X2_List[(N2+1):(N2+N2test)])
TestY<-c(rep(1,N1test),rep(0,N2test))
Testnivector<-rep(ni,(N1test+N2test))

#####################################################
# data preproses
TrainSi_List<-dataSi_f(X_List=TrainX_List,nivector=Trainnivector,p)
TestSi_List<-dataSi_f(X_List=TestX_List,nivector=Testnivector,p)

# DSL training
# Opt_method="RCplex"
Opt_method="Rsymphony"
Constraints="NoConstraints"
# structure='binarystructure'
structure='weightedstructure'
betamethod='majorizationminimization'

h=100;
lamda=1;
mu=1e-4
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

# Test
predprob<-Prediction_f(best_pos,TestSi_List=TestSi_List,Testnivector=Testnivector,
                            AllTheta1=EstimatedTheta1,AllTheta2=EstimatedTheta2,p,h1,h2)
predprob

      
      