######   Solve Discriminant Subgraph Learning (DSL) (Wang, et al., 2021)
#   Input:
#   N1: sample size of class 1
#   N2: sample size of class 2
#   p:  variable number (or the column/row number of network/graph)
#   sub_p: column/row number of the discriminant subgraph
#   h1/h2: hyper-parameters of the Wishart distribution in Class 1/class 2
#   ni: length in each feature matrix xi
#   th: threshold for abs(theta)<th to make Theta to be sparse
#   s: the percentage of non-zero entries picked from subgraph of class 1 to be zero as that of class 1; also the same number of zero entries in the subgraph of class 1 to be nonzero in the subgraph of class 2
   
#   Output:
#   Theta1: shared graph in class 1
#   Theta2: shared graph in class 2
#   sub_pos: positions of subgraph
#   X1_List: data list in Class 1
#   X2_List: data list in Class 2

# Written by Lujia Wang @ Georgia Tech

# library("mvtnorm")
rm(list=ls(all=TRUE));
root="C:\\Users\\lwang724\\Dropbox (GaTech)\\IC code\\Intergration\\github\\"
source(paste(root,"Steps_f.R",sep=''))
source(paste(root,"DataProcedure_f.R",sep=''))
########################3
ni=100
N1=150
N2=N1
p=50;sub_p=5
h1=100
h2=h1
s=0.5
th=0.5


# Theta 1; Theta 2; positions of subgraph
Theta1procedure_result<-Theta1procedure_f(p,th)
Theta1=Theta1procedure_result$Theta1;
Theta2procedure_result<-Theta2procedure_f(Theta1procedure_result,p,sub_p,s,th)
Theta2=Theta2procedure_result$Theta2;
sub_pos=Theta2procedure_result$sub_pos

# Generate the data for each subject in class 1/class 2
X1_result<-Xprocedure_f(Theta=Theta1,p,ni,N=N1,h=h1)
X1_List=X1_result$xList
X2_result<-Xprocedure_f(Theta=Theta2,p,ni,N=N2,h=h2)
X2_List=X2_result$xList
